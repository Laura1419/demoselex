#!/usr/bin/env python3
"""First-pass live reanalysis: our within-trajectory statistic on real AADR genotypes.

For Akbari's tested variants (a time-slope Z / "X" statistic), compute OUR
within-trajectory replication statistic on the real AADR v66.p1 1240K genotypes,
with an empirical null from matched control AADR sites, so the result can be
scattered against Akbari's X.

Design (mirrors docs/I-reanalysis-first-pass.md)
------------------------------------------------
* Test set       : Akbari variants with ``|X| >= 3.16`` (FDR-50% absolute
                   two-sided cut from the paper) intersected with AADR 1240K
                   sites that also appear in the Akbari table.
* Background pool: AADR 1240K sites present in the Akbari table with
                   ``|X| < 3.16`` (absolute two-sided; NOT "X near zero" and
                   NOT the sign of X), randomly sampled with a fixed seed and
                   matched to candidates on chromosome + Akbari allele
                   frequency via a documented caliper matcher (see below).
* Estimator      : ``sim.trajectory_stat.trajectory_statistic`` (sample-size
                   weighted frequency slope + correlated fixed-effect GLS
                   meta-analysis), unchanged.
* Relatedness    : ``data/derived/traj_relatedness.tsv`` is absent, so the scan
                   runs with ``relatedness=None`` (working independence / naive).
                   No relatedness matrix is hand-set.
* Empirical null : a documented caliper matcher, because ``sim/empirical_null.py``
                   requires prespecified metadata features (recombination rate,
                   local LD, imputation quality) that are not available for this
                   AADR panel. Matcher: same chromosome, |pos difference| >= 1 Mb,
                   |AF difference| <= caliper, |missingness difference| <= caliper,
                   |temporal-coverage difference| <= caliper; calipers widen once
                   if fewer than ``--min-controls`` controls survive. Candidates
                   are excluded from every control pool. Per-candidate empirical
                   upper-tail P is ``(1 + #{control score >= candidate score})/(m+1)``.

Negative results are first-class: the script reports what it can and writes
incrementally so a timeout loses nothing.

Usage
-----
    .venv/bin/python analysis/reanalysis_scan.py --max-sites 20000
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "analysis"))
sys.path.insert(0, str(REPO / "sim"))

from load_aadr import (  # noqa: E402
    BUNDLE,
    HEADER_SIZE,
    MAIN_PREFIX,
    GenotypeMatrix,
    parse_header,
)
import load_akbari as ak  # noqa: E402
from trajectory_stat import trajectory_statistic  # noqa: E402

SIG_ABS_X = 3.16  # FDR-50% absolute two-sided cut (paper). Genome-wide cut is 5.45.
BIN_YEARS = 1000
SEED = 42

# Final TSV columns (estimator first; null columns filled in the final rewrite).
ESTIMATOR_COLUMNS = [
    "locus", "source", "chrom", "pos", "akbari_X", "akbari_AF", "missingness",
    "temporal_coverage", "meta_slope", "meta_se", "z_score", "concordance",
    "family_concordance", "n_trajectories", "n_effective", "q", "i_squared",
    "pooled_slope", "is_candidate", "replication_score", "n_controls",
    "p_empirical", "percentile", "bh_fdr",
]


# --------------------------------------------------------------------------- #
# Trajectory cohort
# --------------------------------------------------------------------------- #
def load_cohort(tr_path: Path):
    """Return kept trajectory individuals mapped to AADR matrix rows.

    Columns: ``individual``, ``date_BP``, ``trajectory_id``. Individuals with
    ``trajectory_id < 0`` (unassigned) are dropped.
    """
    tr = pd.read_csv(tr_path, sep="\t",
                     usecols=["individual", "date_BP", "trajectory_id"])
    tr = tr[tr["trajectory_id"] >= 0].copy()
    ind_ids = []
    with open(REPO / BUNDLE / (MAIN_PREFIX + ".ind"), encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if parts:
                ind_ids.append(parts[0])
    id_to_row = {name: i for i, name in enumerate(ind_ids)}
    kept = tr[tr["individual"].isin(id_to_row)].copy()
    n_missing = len(tr) - len(kept)
    kept["row"] = kept["individual"].map(id_to_row).astype(np.int64)
    kept["bin"] = (kept["date_BP"] // BIN_YEARS).astype(np.int64)
    kept["time"] = kept["bin"] * BIN_YEARS + BIN_YEARS // 2  # bin centre, years BP
    kept = kept.sort_values("row").reset_index(drop=True)
    return kept, n_missing


# --------------------------------------------------------------------------- #
# AADR .snp coordinates + Akbari X/AF join
# --------------------------------------------------------------------------- #
def load_snp(path: Path):
    chrom, pos, ref, alt = [], [], [], []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                c = int(parts[1])
            except ValueError:
                continue
            if c < 1 or c > 22:  # Akbari covers autosomes 1-22 only
                continue
            chrom.append(c)
            pos.append(int(parts[3]))
            ref.append(parts[4].upper())
            alt.append(parts[5].upper())
    return (np.asarray(chrom, np.int64), np.asarray(pos, np.int64),
            np.asarray(ref, dtype="U1"), np.asarray(alt, dtype="U1"))


def key(chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    return chrom.astype(np.int64) * 1_000_000_000 + pos.astype(np.int64)


def join_akbari(aadr_key, aadr_ref, aadr_alt, verbose=True):
    """Single streaming pass: fill Akbari X and AF onto each matching AADR site."""
    n = len(aadr_key)
    order = np.argsort(aadr_key, kind="stable")
    sorted_key = aadr_key[order]
    x_arr = np.full(n, np.nan)
    af_arr = np.full(n, np.nan)
    n_rows = 0
    for chunk in ak.iter_chunks():
        n_rows += len(chunk)
        ch = chunk["CHROM"].to_numpy(dtype=np.int64)
        po = chunk["POS"].to_numpy(dtype=np.int64)
        sk = key(ch, po)
        idx = np.searchsorted(sorted_key, sk)
        idx = np.clip(idx, 0, len(sorted_key) - 1)
        hit = sorted_key[idx] == sk
        if not hit.any():
            continue
        aadr_idx = order[idx[hit]]
        ref_ok = aadr_ref[aadr_idx] == chunk["REF"].to_numpy(dtype="U1")[hit]
        alt_ok = aadr_alt[aadr_idx] == chunk["ALT"].to_numpy(dtype="U1")[hit]
        ok = ref_ok & alt_ok
        aadr_idx = aadr_idx[ok]
        x_arr[aadr_idx] = chunk["X"].to_numpy(dtype=float)[hit][ok]
        af_arr[aadr_idx] = chunk["AF"].to_numpy(dtype=float)[hit][ok]
        if verbose and n_rows % 2_000_000 < 1_000_000:
            print(f"  ...akbari rows scanned={n_rows}", flush=True)
    return x_arr, af_arr, n_rows


# --------------------------------------------------------------------------- #
# Genotype panel read (byte-random-access; reads only requested sites)
# --------------------------------------------------------------------------- #
def read_rows_sites(gm, rows, sites):
    """Return uint8 codes (len(rows), len(sites)); 3 = missing, 0/1/2 = alt count."""
    sites = np.asarray(sites, np.int64)
    byte_idx = sites // 4
    shift = ((sites % 4) * 2).astype(np.uint8)
    out = np.empty((len(rows), len(sites)), dtype=np.uint8)
    mm = gm._mm  # noqa: SLF001 - the reader exposes no byte-sliced API
    for k, i in enumerate(rows):
        out[k] = (mm[int(i), byte_idx] >> shift) & np.uint8(0x3)
    return out


# --------------------------------------------------------------------------- #
# Main scan
# --------------------------------------------------------------------------- #
def build_panel(max_sites, verbose=True):
    """Choose candidate + background AADR site indices within the site cap."""
    aadr_dir = REPO / BUNDLE
    snp_path = aadr_dir / (MAIN_PREFIX + ".snp")
    if verbose:
        print("loading AADR .snp coordinates ...", flush=True)
    chrom, pos, ref, alt = load_snp(snp_path)
    akey = key(chrom, pos)
    if verbose:
        print(f"  AADR autosome sites (1-22): {len(akey)}", flush=True)
        print("joining Akbari X/AF (streaming 9.7M rows) ...", flush=True)
    x_arr, af_arr, n_ak_rows = join_akbari(akey, ref, alt, verbose=verbose)
    have = np.isfinite(x_arr)
    sig_mask = have & (np.abs(x_arr) >= SIG_ABS_X)
    bg_mask = have & (np.abs(x_arr) < SIG_ABS_X)

    sig_idx = np.nonzero(sig_mask)[0]
    bg_idx = np.nonzero(bg_mask)[0]
    rng = np.random.default_rng(SEED)
    n_cand = min(len(sig_idx), max(1000, max_sites // 4))
    if len(sig_idx) > n_cand:
        sig_idx = np.sort(rng.choice(sig_idx, size=n_cand, replace=False))
    # Shuffle background so non-significant rows are an unbiased pool, then downsample.
    bg_idx = bg_idx[rng.permutation(len(bg_idx))]
    n_bg = max(0, max_sites - len(sig_idx))
    bg_idx = np.sort(bg_idx[:n_bg])

    panel = np.concatenate([sig_idx, bg_idx])
    is_cand = np.zeros(len(panel), dtype=bool)
    is_cand[:len(sig_idx)] = True
    stats = {
        "n_aadr_autosome_sites": int(len(akey)),
        "n_akbari_rows_scanned": int(n_ak_rows),
        "n_aadr_sites_in_akbari": int(have.sum()),
        "n_akbari_significant_absX_ge_3p16_in_aadr": int(len(np.nonzero(sig_mask)[0])),
        "n_background_pool_in_aadr": int(len(np.nonzero(bg_mask)[0])),
        "n_candidates_analysed": int(len(sig_idx)),
        "n_background_analysed": int(len(bg_idx)),
        "n_sites_analysed": int(len(panel)),
    }
    labels = np.array([f"{chrom[i]}:{pos[i]}:{ref[i]}:{alt[i]}" for i in panel])
    return {
        "panel": panel, "is_cand": is_cand, "labels": labels,
        "chrom": chrom[panel], "pos": pos[panel],
        "akbari_X": x_arr[panel], "akbari_AF": af_arr[panel], "stats": stats,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--max-sites", type=int, default=30000,
                    help="total cap on candidate+background sites analysed")
    ap.add_argument("--chunk-sites", type=int, default=2000)
    ap.add_argument("--min-controls", type=int, default=20)
    ap.add_argument("--af-caliper", type=float, default=0.10)
    ap.add_argument("--miss-caliper", type=float, default=0.10)
    ap.add_argument("--cov-caliper", type=float, default=0.15)
    ap.add_argument("--out-tsv", type=Path,
                    default=REPO / "data/derived/replication_stats.tsv")
    ap.add_argument("--out-json", type=Path,
                    default=REPO / "data/derived/reanalysis_metrics.json")
    args = ap.parse_args(argv)

    t0 = time.time()
    relatedness_path = REPO / "data/derived/traj_relatedness.tsv"
    rel = None
    relatedness_used = relatedness_path.exists()
    if relatedness_used:
        rel = pd.read_csv(relatedness_path, sep="\t", index_col=0)
        rel.index = [f"t{i}" for i in rel.index]
        rel.columns = [f"t{i}" for i in rel.columns]
        print(f"relatedness: applying {relatedness_path.name} ({rel.shape[0]}x"
              f"{rel.shape[1]}) to trajectory_statistic.", flush=True)
    else:
        print("WARNING: data/derived/traj_relatedness.tsv absent -> relatedness=None "
              "(naive working independence). No matrix hand-set.", flush=True)

    cohort, n_unmapped = load_cohort(REPO / "data/derived/trajectories.tsv")
    rows = np.unique(cohort["row"].to_numpy())
    row_pos = {int(r): k for k, r in enumerate(rows)}
    cohort_pos = cohort["row"].map(row_pos).to_numpy(np.int64)
    # Map each kept individual to a (trajectory, bin) group.
    gkeys = (cohort["trajectory_id"].to_numpy() * 100_000
             + cohort["bin"].to_numpy())
    uniq_g, g_of_row = np.unique(gkeys, return_inverse=True)
    n_groups = len(uniq_g)
    group_traj = (uniq_g // 100_000).astype(np.int64)
    group_bin = (uniq_g % 100_000).astype(np.int64)
    group_time = group_bin * BIN_YEARS + BIN_YEARS // 2
    row_group = np.empty(len(rows), np.int64)
    row_group[cohort_pos] = g_of_row
    sorter = np.argsort(row_group, kind="stable")
    gr_sorted = row_group[sorter]
    starts = np.searchsorted(gr_sorted, np.arange(n_groups))

    print(f"cohort: {len(cohort)} individuals, {len(rows)} AADR rows, "
          f"{n_groups} (trajectory,bin) cells; unmapped={n_unmapped}", flush=True)

    panel_info = build_panel(args.max_sites)
    panel = panel_info["panel"]
    is_cand = panel_info["is_cand"]
    labels = panel_info["labels"]
    n_panel = len(panel)
    n_ind = len(rows)

    gm = GenotypeMatrix(str(REPO / BUNDLE / (MAIN_PREFIX + ".geno")))
    assert gm.n_ind == parse_header(gm.path).n_ind

    # Incremental checkpoint writer (estimator columns only until the final pass).
    out_tsv = args.out_tsv
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", encoding="utf-8") as fh:
        fh.write("\t".join(ESTIMATOR_COLUMNS) + "\n")

    summaries = []
    miss = np.full(n_panel, np.nan)
    tcov = np.full(n_panel, np.nan)
    n_dropped_support = 0
    for start in range(0, n_panel, args.chunk_sites):
        stop = min(start + args.chunk_sites, n_panel)
        sites = panel[start:stop]
        codes = read_rows_sites(gm, rows, sites)  # (n_ind, n_site_chunk)
        alt_mat = np.where(codes <= 2, codes, 0).astype(np.int32)
        tot_mat = (2 * (codes <= 2)).astype(np.int32)
        alt_sum = np.add.reduceat(alt_mat[sorter], starts, axis=0)  # (n_groups, n_chunk)
        tot_sum = np.add.reduceat(tot_mat[sorter], starts, axis=0)
        miss[start:stop] = 1.0 - tot_sum.sum(axis=0) / (2.0 * n_ind)
        tcov[start:stop] = (tot_sum > 0).sum(axis=0) / float(n_groups)

        gi, si = np.nonzero(tot_sum)
        if len(gi) == 0:
            n_dropped_support += stop - start
            continue
        tidy = pd.DataFrame({
            "locus": labels[start + si],
            "trajectory": np.char.add("t", group_traj[gi].astype(str)),
            "time": group_time[gi].astype(float),
            "alt_count": alt_sum[gi, si].astype(float),
            "total_count": tot_sum[gi, si].astype(float),
        })
        # Keep only (locus, trajectory) with >= 2 distinct times, then loci with
        # >= 2 informative trajectories (required by the estimator).
        pair = tidy["locus"] + "\x1f" + tidy["trajectory"]
        per = tidy.assign(pair=pair).groupby("pair", sort=False)["time"].nunique()
        keep_pairs = set(per[per >= 2].index)
        tidy = tidy[pair.isin(keep_pairs)]
        ntraj = tidy.groupby("locus", sort=False)["trajectory"].nunique()
        good_loci = set(ntraj[ntraj >= 2].index)
        dropped = stop - start - len(good_loci)
        n_dropped_support += dropped
        tidy = tidy[tidy["locus"].isin(good_loci)]
        if tidy.empty:
            continue
        res = trajectory_statistic(tidy, relatedness=rel)
        summ = res.summary
        chunk_labels = labels[start:stop]
        chunk_cand = is_cand[start:stop]
        meta = pd.DataFrame({
            "locus": chunk_labels, "chrom": panel_info["chrom"][start:stop],
            "pos": panel_info["pos"][start:stop],
            "akbari_X": panel_info["akbari_X"][start:stop],
            "akbari_AF": panel_info["akbari_AF"][start:stop],
            "missingness": miss[start:stop],
            "temporal_coverage": tcov[start:stop],
            "is_candidate": chunk_cand})
        summ = summ.merge(meta, on="locus", how="left")
        summaries.append(summ)

        chunk_out = summ.copy()
        chunk_out["source"] = np.where(chunk_out["is_candidate"],
                                       "akbari_sig", "background")
        chunk_out["replication_score"] = (chunk_out["z_score"].abs()
                                          * chunk_out["family_concordance"])
        for col in ESTIMATOR_COLUMNS:
            if col not in chunk_out:
                chunk_out[col] = np.nan
        with open(out_tsv, "a", encoding="utf-8") as fh:
            chunk_out[ESTIMATOR_COLUMNS].to_csv(fh, sep="\t", index=False,
                                                header=False)
        print(f"  chunk {start}-{stop}: {len(summ)} loci scored, "
              f"{dropped} dropped (<2 informative trajectories); t={time.time()-t0:.0f}s",
              flush=True)

    summary = pd.concat(summaries, ignore_index=True) if summaries else pd.DataFrame()
    n_scored = len(summary)
    # Assertion required by the task validation.
    assert (summary["n_trajectories"] >= 2).all(), "a scored locus had < 2 trajectories"
    print(f"scored loci={n_scored}; dropped for <2 informative trajectories="
          f"{n_dropped_support}", flush=True)

    # ------------------------------------------------------------------ null #
    summary["replication_score"] = (summary["z_score"].abs()
                                    * summary["family_concordance"])
    finite = np.isfinite(summary["replication_score"].to_numpy(dtype=float))
    cand = summary[summary["is_candidate"] & finite].copy()
    bg = summary[(~summary["is_candidate"]) & finite].copy()

    p_emp = np.full(len(cand), np.nan)
    n_ctrl = np.zeros(len(cand), dtype=int)
    pct = np.full(len(cand), np.nan)
    exact = np.zeros(len(cand), dtype=bool)
    for j, (cidx, crow) in enumerate(cand.iterrows()):
        pool = bg[(bg["chrom"] == crow["chrom"])
                  & ((bg["pos"] - crow["pos"]).abs() >= 1_000_000)]
        chosen = None
        for scale in (1.0, 2.0, 4.0):
            af_c = min(args.af_caliper * scale, 0.5)
            mi_c = min(args.miss_caliper * scale, 0.5)
            co_c = min(args.cov_caliper * scale, 0.5)
            sel = pool[((pool["akbari_AF"] - crow["akbari_AF"]).abs() <= af_c)
                       & ((pool["missingness"] - crow["missingness"]).abs() <= mi_c)
                       & ((pool["temporal_coverage"]
                           - crow["temporal_coverage"]).abs() <= co_c)]
            if len(sel) >= args.min_controls:
                chosen = sel
                exact[j] = scale == 1.0
                break
        if chosen is None:
            chosen = pool  # documented last-resort fallback: all same-chrom controls
        m = len(chosen)
        scores = chosen["replication_score"].to_numpy(dtype=float)
        cs = float(crow["replication_score"])
        n_ctrl[j] = m
        if m == 0:
            continue
        less = int((scores < cs).sum())
        equal = int((scores == cs).sum())
        p_emp[j] = (1 + m - less) / (m + 1)
        pct[j] = 100.0 * (less + 0.5 * equal) / m

    cand["p_empirical"] = p_emp
    cand["n_controls"] = n_ctrl
    cand["percentile"] = pct
    cand["matched_exact"] = exact

    # Benjamini-Hochberg FDR over the Akbari test set (finite p only).
    cand["bh_fdr"] = np.nan
    ok = np.isfinite(cand["p_empirical"].to_numpy(dtype=float))
    if ok.any():
        pv = cand.loc[ok, "p_empirical"].to_numpy(dtype=float)
        order = np.argsort(pv)
        ranked = pv[order]
        m = len(ranked)
        q = ranked * m / np.arange(1, m + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)
        out = np.empty(m)
        out[order] = q
        cand.loc[ok, "bh_fdr"] = out

    # Assemble the full table: every scored locus is kept; the null columns are
    # filled for candidate loci (background rows keep NaN there).
    summary = summary.copy()
    for col in ("replication_score", "n_controls", "p_empirical", "percentile",
                "bh_fdr"):
        if col not in summary:
            summary[col] = np.nan
    summary.loc[cand.index, "n_controls"] = n_ctrl
    summary.loc[cand.index, "p_empirical"] = p_emp
    summary.loc[cand.index, "percentile"] = pct
    summary.loc[cand.index, "bh_fdr"] = cand["bh_fdr"].to_numpy()
    if "source" not in summary:
        summary["source"] = "scan"
    summary = summary[ESTIMATOR_COLUMNS]
    summary.to_csv(out_tsv, sep="\t", index=False)

    # ---------------------------------------------------------------- output #
    sig_scores = summary.loc[summary["is_candidate"] & finite, "replication_score"]
    bg_scores = summary.loc[~summary["is_candidate"] & finite, "replication_score"]
    metrics = dict(panel_info["stats"], **{
        "seed": SEED,
        "bin_years": BIN_YEARS,
        "sig_abs_x_threshold": SIG_ABS_X,
        "genome_wide_abs_x": 5.45,
        "n_trajectory_individuals": int(len(cohort)),
        "n_aadr_rows_used": int(n_ind),
        "n_trajectory_bin_cells": int(n_groups),
        "n_loci_scored": int(n_scored),
        "n_dropped_lt2_informative_trajectories": int(n_dropped_support),
        "n_candidates_scored": int(summary["is_candidate"].sum()),
        "n_background_scored": int((~summary["is_candidate"]).sum()),
        "n_background_scored_finite_score": int(len(bg)),
        "relatedness_applied": bool(relatedness_used),
        "null_method": "caliper matcher (chrom, AF, missingness, temporal coverage; "
                       "same-chrom >=1Mb excluded; calipers widen 1x/2x/4x)",
        "min_controls": args.min_controls,
        "n_candidates_matched_exactly": int(cand["matched_exact"].sum()),
        "n_candidates_finite_p": int(np.isfinite(cand["p_empirical"]).sum()),
        "n_candidates_bh_fdr_le_0p05": int((cand["bh_fdr"] <= 0.05).sum()),
        "candidate_replication_score_mean": float(sig_scores.mean()) if len(sig_scores) else None,
        "background_replication_score_mean": float(bg_scores.mean()) if len(bg_scores) else None,
        "candidate_replication_score_median": float(sig_scores.median()) if len(sig_scores) else None,
        "background_replication_score_median": float(bg_scores.median()) if len(bg_scores) else None,
        "runtime_seconds": round(time.time() - t0, 1),
    })
    args.out_json.write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"wrote {out_tsv} ({len(summary)} rows) and {args.out_json}", flush=True)
    print(json.dumps({k: metrics[k] for k in (
        "n_sites_analysed", "n_candidates_scored", "n_background_scored",
        "n_loci_scored", "n_dropped_lt2_informative_trajectories",
        "relatedness_applied", "runtime_seconds")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
