#!/usr/bin/env python3
"""Tier-2 decisive like-for-like test: OUR trajectory statistic on Akbari's own
imputed genotypes, same individuals and same variants as Akbari et al. 2026.

This replaces the Tier-1 AADR 1240K subset (`analysis/reanalysis_scan.py`) with
the Akbari PLINK2 imputed bundle, keeping the individual set and variant set the
paper actually analysed.

Design
------
* Samples  : the pgen indices of the trajectory individuals present in the bundle
             (explicit id -> pgen-index mapping; report the count).
* Variants : panel = (Akbari |X| >= 3.16) intersect the HighQualityVariants QC
             list, plus a background pool (QC, |X| < 3.16) used as matched
             controls. The panel is capped to fit the budget; both the cap and
             the pre-cap availability are reported.
* Units    : the six merged replication groups of `data/derived/trajectories_v3.tsv`
             (`merged_group` G0..G5), with the matching `data/derived/traj_relatedness_v3.tsv`
             (6x6) matrix. Units and matrix match; no matrix is hand-set. See
             docs/traj-v3-merge.md: the merged groups are the genetically-independent
             replication units (13 v2 trajectories are over-split). The G0 bloc is
             known to be an over-merge -> reported as a limitation.
* Statistics: BOTH
    (a) the standard replication statistic via
        `sim.trajectory_stat.trajectory_statistic` (relatedness applied);
    (b) the PRE-REGISTERED max-over-subsets statistic (max |z| over subsets of
        size 2..4 of the present groups), reusing the GLS z of
        `analysis/subset_regime.py` (docs/subset-regime-prereg.md).
* Null     : matched-control empirical null for EACH statistic separately, using
             the caliper matcher of `analysis/reanalysis_scan.py` (chrom, Akbari
             AF, missingness, temporal coverage; same-chrom >= 1 Mb excluded;
             calipers widen 1x/2x/4x). BH-FDR across tested variants.

Writes `data/derived/tier2_stats.tsv` incrementally after every chunk (header
first) so a timeout loses nothing. Negative results are first-class.

Usage
-----
    .venv/bin/python analysis/tier2_scan.py --max-variants 2000
"""

from __future__ import annotations

import argparse
import itertools
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "analysis"))
sys.path.insert(0, str(REPO / "sim"))

import load_akbari as ak  # noqa: E402
import load_akbari_geno as akg  # noqa: E402
import reanalysis_scan as rs  # noqa: E402
from subset_regime import gls_z  # noqa: E402  (reuse the pre-registered GLS z)
from trajectory_stat import _correlation, _fit_slopes, trajectory_statistic  # noqa: E402

SIG_ABS_X = 3.16
GW_ABS_X = 5.45
BIN_YEARS = 1000
SEED = 42
MIN_SUBSET = 2
MAX_SUBSET = 4

QC_PATH = REPO / "data/raw/akbari/tier2/HighQualityVariants.tsv"
TRAJ_V3 = REPO / "data/derived/trajectories_v3.tsv"
REL_V3 = REPO / "data/derived/traj_relatedness_v3.tsv"
POOL_CACHE = Path("/tmp/tier2_pool_cache.pkl")
MAP_CACHE = Path("/tmp/tier2_sample_map.pkl")

# TSV columns. The estimator block is written incrementally; the null columns for
# the two statistics are filled in the final rewrite.
COLUMNS = [
    "locus", "chrom", "pos", "akbari_X", "akbari_AF", "missingness",
    "temporal_coverage", "n_groups", "meta_slope", "meta_se", "z_score",
    "concordance", "family_concordance", "n_effective", "q", "i_squared",
    "pooled_slope", "subset_max_stat", "argmax_subset_size", "argmax_subset",
    "argmax_z", "complement_slope", "complement_z", "complement_flat",
    "complement_opposite", "n_subsets", "is_candidate", "source",
    "replication_score", "subset_score", "n_controls", "p_empirical", "bh_fdr",
    "n_controls_subset", "p_empirical_subset", "bh_fdr_subset",
]


# --------------------------------------------------------------------------- #
# Sample mapping: trajectory individuals -> pgen sample indices
# --------------------------------------------------------------------------- #
def build_sample_map(verbose=True):
    """Return a dataframe of kept individuals with columns ``individual``,
    ``merged_group``, ``bin``, ``time``, ``pgen_index`` (position within the
    requested sample subset). Unassigned individuals are dropped."""
    if MAP_CACHE.exists():
        return pickle.loads(MAP_CACHE.read_bytes())
    sample_ids = akg.parse_psam(akg.DEFAULT_DIR / f"{akg.PREFIX}.psam")
    id_to_full = {s: i for i, s in enumerate(sample_ids)}
    tr = pd.read_csv(TRAJ_V3, sep="\t",
                     usecols=["individual", "date_BP", "trajectory_id", "merged_group"])
    tr = tr[tr["merged_group"] != "unassigned"].copy()
    tr = tr[tr["individual"].isin(id_to_full)].copy()
    tr["pgen_full_index"] = tr["individual"].map(id_to_full).astype(np.int64)
    tr = tr.sort_values("pgen_full_index").reset_index(drop=True)
    tr["pgen_index"] = np.arange(len(tr), dtype=np.int64)  # position in subset vec
    tr["bin"] = (tr["date_BP"] // BIN_YEARS).astype(np.int64)
    tr["time"] = tr["bin"] * BIN_YEARS + BIN_YEARS // 2
    if verbose:
        print(f"sample map: {len(tr)} trajectory individuals present in pgen; "
              f"groups={tr['merged_group'].value_counts().to_dict()}", flush=True)
    MAP_CACHE.write_bytes(pickle.dumps(tr))
    return tr


# --------------------------------------------------------------------------- #
# Panel pool: Akbari |X| bands intersect the QC list, with pgen variant indices
# --------------------------------------------------------------------------- #
def build_pool(max_bg, verbose=True):
    """Return a dataframe of QC-confirmed Akbari variants with pgen indices.

    Columns: ID, CHROM, POS, REF, ALT, AF, X, is_candidate, pgen_index.
    Candidates = |X| >= 3.16; background = |X| < 3.16 (sampled down to max_bg).
    ``pgen_index`` is -1 when the variant is absent from the pgen bundle.
    """
    if POOL_CACHE.exists():
        return pickle.loads(POOL_CACHE.read_bytes())

    t0 = time.time()
    if verbose:
        print("loading HighQualityVariants QC id set ...", flush=True)
    qc_ids = set()
    with open(QC_PATH, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                qc_ids.add(fields[2])  # ID column: CHROM POS ID REF ALT
    if verbose:
        print(f"  QC variants: {len(qc_ids):,} ({time.time()-t0:.0f}s)", flush=True)

    rng = np.random.default_rng(SEED)
    cand_rows, bg_rows = [], []
    n_rows = 0
    for chunk in ak.iter_chunks(usecols=["ID", "CHROM", "POS", "REF", "ALT", "AF", "X"]):
        n_rows += len(chunk)
        chunk = chunk[chunk["ID"].isin(qc_ids)]
        if not len(chunk):
            continue
        absx = chunk["X"].abs().to_numpy()
        c = chunk[absx >= SIG_ABS_X]
        if len(c):
            cand_rows.append(c)
        b = chunk[absx < SIG_ABS_X]
        if len(b):
            bg_rows.append(b)
    cand = (pd.concat(cand_rows, ignore_index=True) if cand_rows
            else pd.DataFrame(columns=["ID", "CHROM", "POS", "REF", "ALT", "AF", "X"]))
    bg = (pd.concat(bg_rows, ignore_index=True) if bg_rows
          else pd.DataFrame(columns=cand.columns))
    n_cand_avail, n_bg_avail = len(cand), len(bg)
    # Sample the background pool so the pool build stays bounded.
    if len(bg) > max_bg:
        bg = bg.iloc[np.sort(rng.choice(len(bg), size=max_bg, replace=False))].copy()
    if verbose:
        print(f"  Akbari rows scanned={n_rows:,}; candidates(|X|>={SIG_ABS_X})={n_cand_avail:,}; "
              f"background(|X|<{SIG_ABS_X})={n_bg_avail:,} -> pool {len(bg):,} "
              f"({time.time()-t0:.0f}s)", flush=True)

    pool = pd.concat([cand, bg], ignore_index=True)
    pool["is_candidate"] = np.r_[np.ones(len(cand), bool), np.zeros(len(bg), bool)]
    pool = pool.drop_duplicates("ID").reset_index(drop=True)

    # Map each pool variant ID to its pgen variant index (single streaming pass).
    if verbose:
        print("streaming pvar to map panel ids -> pgen variant indices ...", flush=True)
    needed = set(pool["ID"])
    id_to_index = {}
    for i, (_chrom, _pos, vid, _ref, _alt) in enumerate(akg.iter_pvar(akg.DEFAULT_DIR / f"{akg.PREFIX}.pvar")):
        if vid in needed:
            id_to_index[vid] = i
    pool["pgen_index"] = pool["ID"].map(id_to_index).fillna(-1).astype(np.int64)
    n_missing = int((pool["pgen_index"] < 0).sum())
    if verbose:
        print(f"  mapped {len(id_to_index):,}/{len(pool):,} pool ids; absent={n_missing:,} "
              f"({time.time()-t0:.0f}s)", flush=True)
    pool.attrs["n_cand_avail"] = n_cand_avail
    pool.attrs["n_bg_avail"] = n_bg_avail
    pool.attrs["n_akbari_rows_scanned"] = n_rows
    pool.attrs["n_qc"] = len(qc_ids)
    POOL_CACHE.write_bytes(pickle.dumps(pool))
    return pool


def select_panel(pool, max_variants):
    """Cap the pool: up to max_variants//4 candidates and the rest background."""
    pool = pool[pool["pgen_index"] >= 0].reset_index(drop=True)
    cand_idx = np.nonzero(pool["is_candidate"].to_numpy())[0]
    bg_idx = np.nonzero(~pool["is_candidate"].to_numpy())[0]
    rng = np.random.default_rng(SEED)
    n_cand = min(len(cand_idx), max(1000, max_variants // 4))
    if len(cand_idx) > n_cand:
        cand_idx = np.sort(rng.choice(cand_idx, size=n_cand, replace=False))
    n_bg = min(len(bg_idx), max(0, max_variants - len(cand_idx)))
    if len(bg_idx) > n_bg:
        bg_idx = np.sort(rng.choice(bg_idx, size=n_bg, replace=False))
    take = np.concatenate([cand_idx, bg_idx])
    panel = pool.iloc[take].reset_index(drop=True)
    panel = panel.sort_values("pgen_index").reset_index(drop=True)
    return panel


# --------------------------------------------------------------------------- #
# Pre-registered max-over-subsets statistic (reuses subset_regime.gls_z)
# --------------------------------------------------------------------------- #
def merged_subset_statistic(slopes_locus, r_full, group_order):
    """Max |z| over subsets of size 2..4 of the present groups (pre-registered).

    ``slopes_locus`` maps group label -> (slope, se); ``group_order`` is the fixed
    group label ordering used to index ``r_full``.
    """
    present = [g for g in group_order if g in slopes_locus]
    gidx = {g: i for i, g in enumerate(group_order)}
    best = None
    n_sub = 0
    for size in range(MIN_SUBSET, MAX_SUBSET + 1):
        for combo in itertools.combinations(present, size):
            n_sub += 1
            ridx = [gidx[g] for g in combo]
            r_sub = r_full[np.ix_(ridx, ridx)]
            b = np.array([slopes_locus[g][0] for g in combo], float)
            s = np.array([slopes_locus[g][1] for g in combo], float)
            meta, z = gls_z(b, s, r_sub)
            if best is None or abs(z) > abs(best[2]):
                best = (combo, meta, z)
    if best is None:
        return None
    combo, meta, z = best
    comp = [g for g in present if g not in combo]
    comp_slope = comp_z = np.nan
    comp_flat = comp_opp = False
    if len(comp) >= 2:
        ridx = [gidx[g] for g in comp]
        r_sub = r_full[np.ix_(ridx, ridx)]
        b = np.array([slopes_locus[g][0] for g in comp], float)
        s = np.array([slopes_locus[g][1] for g in comp], float)
        c_meta, c_z = gls_z(b, s, r_sub)
        comp_slope, comp_z = c_meta, c_z
        comp_flat = abs(c_z) < 1.0
        comp_opp = np.sign(c_meta) != np.sign(meta) and c_meta != 0 and meta != 0
    return {
        "n_groups": len(present), "n_subsets": n_sub,
        "subset_max_stat": abs(z), "argmax_subset_size": len(combo),
        "argmax_subset": ",".join(combo), "argmax_z": z,
        "complement_slope": comp_slope, "complement_z": comp_z,
        "complement_flat": bool(comp_flat), "complement_opposite": bool(comp_opp),
    }


# --------------------------------------------------------------------------- #
# Matched-control empirical null + BH-FDR (same matcher as reanalysis_scan)
# --------------------------------------------------------------------------- #
def empirical_null(summary, score_col, args):
    finite = np.isfinite(summary[score_col].to_numpy(dtype=float))
    cand = summary[summary["is_candidate"] & finite].copy()
    bg = summary[(~summary["is_candidate"]) & finite].copy()
    p_emp = np.full(len(cand), np.nan)
    n_ctrl = np.zeros(len(cand), dtype=int)
    for j, (_i, crow) in enumerate(cand.iterrows()):
        pool = bg[(bg["chrom"] == crow["chrom"])
                  & ((bg["pos"] - crow["pos"]).abs() >= 1_000_000)]
        chosen = None
        for scale in (1.0, 2.0, 4.0):
            sel = pool[((pool["akbari_AF"] - crow["akbari_AF"]).abs() <= min(args.af_caliper * scale, 0.5))
                       & ((pool["missingness"] - crow["missingness"]).abs() <= min(args.miss_caliper * scale, 0.5))
                       & ((pool["temporal_coverage"] - crow["temporal_coverage"]).abs() <= min(args.cov_caliper * scale, 0.5))]
            if len(sel) >= args.min_controls:
                chosen = sel
                break
        if chosen is None:
            chosen = pool  # documented last-resort fallback
        m = len(chosen)
        n_ctrl[j] = m
        if m == 0:
            continue
        scores = chosen[score_col].to_numpy(dtype=float)
        cs = float(crow[score_col])
        p_emp[j] = (1 + int((scores >= cs).sum())) / (m + 1)
    bh = np.full(len(cand), np.nan)
    ok = np.isfinite(p_emp)
    if ok.any():
        pv = p_emp[ok]
        order = np.argsort(pv)
        ranked = pv[order]
        m = len(ranked)
        q = ranked * m / np.arange(1, m + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q = np.clip(q, 0, 1)
        out = np.empty(m)
        out[order] = q
        bh[ok] = out
    return cand.index, p_emp, n_ctrl, bh


# --------------------------------------------------------------------------- #
# Main scan
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--max-variants", type=int, default=20000,
                    help="total cap on candidate+background variants analysed")
    ap.add_argument("--chunk-variants", type=int, default=2000)
    ap.add_argument("--max-bg-pool", type=int, default=150000)
    ap.add_argument("--min-controls", type=int, default=20)
    ap.add_argument("--af-caliper", type=float, default=0.10)
    ap.add_argument("--miss-caliper", type=float, default=0.10)
    ap.add_argument("--cov-caliper", type=float, default=0.15)
    ap.add_argument("--out-tsv", type=Path,
                    default=REPO / "data/derived/tier2_stats.tsv")
    ap.add_argument("--out-json", type=Path,
                    default=REPO / "data/derived/tier2_metrics.json")
    args = ap.parse_args(argv)

    t0 = time.time()
    smap = build_sample_map()
    group_order = list(pd.read_csv(REL_V3, sep="\t", index_col=0).index)
    rel = pd.read_csv(REL_V3, sep="\t", index_col=0)
    rel.index = group_order
    rel.columns = group_order
    r_full = _correlation(rel, group_order)
    print(f"units: {len(group_order)} merged groups {group_order}; "
          f"relatedness={REL_V3.name} ({rel.shape[0]}x{rel.shape[1]}); "
          f"N_eff={len(group_order)**2 / r_full.sum():.3f}", flush=True)

    pool = build_pool(args.max_bg_pool)
    panel = select_panel(pool, args.max_variants)
    n_panel = len(panel)
    print(f"panel: {n_panel} variants analysed "
          f"(candidates={int(panel['is_candidate'].sum())}, "
          f"background={int((~panel['is_candidate']).sum())}); "
          f"available before cap: candidates={pool.attrs['n_cand_avail']:,}, "
          f"background={pool.attrs['n_bg_avail']:,}", flush=True)

    # Sample vector = pgen_full_index of the kept individuals, in subset order.
    subset = smap["pgen_full_index"].to_numpy(np.int32)
    n_ind = len(subset)
    # (group, bin) cells for the 5788 individuals.
    gord = {g: i for i, g in enumerate(group_order)}
    gid = smap["merged_group"].map(gord).to_numpy(np.int64)
    binv = smap["bin"].to_numpy(np.int64)
    gkeys = gid * 100_000 + binv
    uniq_g, g_of_row = np.unique(gkeys, return_inverse=True)
    n_groups_all = len(uniq_g)
    group_traj = (uniq_g // 100_000).astype(np.int64)
    group_bin = (uniq_g % 100_000).astype(np.int64)
    group_time = group_bin * BIN_YEARS + BIN_YEARS // 2
    sorter = np.argsort(g_of_row, kind="stable")
    gr_sorted = g_of_row[sorter]
    starts = np.searchsorted(gr_sorted, np.arange(n_groups_all))
    labels = (panel["CHROM"].astype(str) + ":" + panel["POS"].astype(str) + ":"
              + panel["REF"].astype(str) + ":" + panel["ALT"].astype(str))

    out_tsv = args.out_tsv
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", encoding="utf-8") as fh:
        fh.write("\t".join(COLUMNS) + "\n")

    geno = akg.AkbariGeno(sample_subset=subset)
    assert geno.n_samples == n_ind

    summaries = []
    miss = np.full(n_panel, np.nan)
    tcov = np.full(n_panel, np.nan)
    n_dropped = 0
    for start in range(0, n_panel, args.chunk_variants):
        stop = min(start + args.chunk_variants, n_panel)
        pgen_idx = panel["pgen_index"].to_numpy()[start:stop]
        dos = geno.read_variants(pgen_idx)  # (chunk, n_ind) float32, NaN missing
        finite = np.isfinite(dos)
        alt_mat = np.where(finite, np.rint(dos), 0).astype(np.int32)
        tot_mat = (2 * finite).astype(np.int32)
        alt_sum = np.add.reduceat(alt_mat[sorter], starts, axis=0)  # (n_groups_all, chunk)
        tot_sum = np.add.reduceat(tot_mat[sorter], starts, axis=0)
        miss[start:stop] = 1.0 - tot_sum.sum(axis=0) / (2.0 * n_ind)
        tcov[start:stop] = (tot_sum > 0).sum(axis=0) / float(n_groups_all)

        chunk_labels = labels[start:stop].to_numpy()
        gi, si = np.nonzero(tot_sum)
        chunk_rows = []
        if len(gi):
            tidy = pd.DataFrame({
                "locus": chunk_labels[si],
                "trajectory": np.array([f"G{group_traj[g]}" for g in gi]),
                "time": group_time[gi].astype(float),
                "alt_count": alt_sum[gi, si].astype(np.int64),
                "total_count": tot_sum[gi, si].astype(np.int64),
            })
            # >= 2 distinct times per (locus, group), then >= 2 informative groups.
            pair = tidy["locus"] + "\x1f" + tidy["trajectory"]
            per = tidy.assign(pair=pair).groupby("pair", sort=False)["time"].nunique()
            tidy = tidy[pair.isin(set(per[per >= 2].index))]
            ntraj = tidy.groupby("locus", sort=False)["trajectory"].nunique()
            good = set(ntraj[ntraj >= 2].index)
            tidy = tidy[tidy["locus"].isin(good)]
            if len(tidy):
                res = trajectory_statistic(tidy, relatedness=rel)
                summ = res.summary.set_index("locus")
                slopes = _fit_slopes(tidy)
                sl_by_locus: dict[str, dict[str, tuple[float, float]]] = {}
                for row in slopes.itertuples(index=False):
                    sl_by_locus.setdefault(row.locus, {})[row.trajectory] = (row.slope, row.se)
                for lab in summ.index:
                    meta = summ.loc[lab]
                    st = merged_subset_statistic(sl_by_locus.get(lab, {}),
                                                 r_full, group_order)
                    chunk_rows.append({
                        "locus": lab,
                        "n_groups": int(meta["n_trajectories"]),
                        "meta_slope": meta["meta_slope"], "meta_se": meta["meta_se"],
                        "z_score": meta["z_score"], "concordance": meta["concordance"],
                        "family_concordance": meta["family_concordance"],
                        "n_effective": meta["n_effective"], "q": meta["q"],
                        "i_squared": meta["i_squared"], "pooled_slope": meta["pooled_slope"],
                        **(st or {k: np.nan for k in (
                            "subset_max_stat", "argmax_subset_size", "argmax_subset",
                            "argmax_z", "complement_slope", "complement_z",
                            "complement_flat", "complement_opposite", "n_subsets")}),
                    })
        scored = pd.DataFrame(chunk_rows)
        n_scored_chunk = len(scored)
        dropped_chunk = (stop - start) - n_scored_chunk
        n_dropped += dropped_chunk
        # Emit every chunk locus; unscored rows keep NaN estimator columns.
        chunk_out = pd.DataFrame({
            "locus": chunk_labels,
            "chrom": panel["CHROM"].to_numpy()[start:stop],
            "pos": panel["POS"].to_numpy()[start:stop],
            "akbari_X": panel["X"].to_numpy()[start:stop],
            "akbari_AF": panel["AF"].to_numpy()[start:stop],
            "missingness": miss[start:stop],
            "temporal_coverage": tcov[start:stop],
            "is_candidate": panel["is_candidate"].to_numpy()[start:stop],
        })
        if len(scored):
            chunk_out = chunk_out.merge(scored, on="locus", how="left")
        for col in COLUMNS:
            if col not in chunk_out:
                chunk_out[col] = np.nan
        chunk_out["source"] = np.where(chunk_out["is_candidate"], "akbari_sig", "background")
        chunk_out["replication_score"] = chunk_out["z_score"].abs() * chunk_out["family_concordance"]
        chunk_out["subset_score"] = chunk_out["subset_max_stat"]
        with open(out_tsv, "a", encoding="utf-8") as fh:
            chunk_out[COLUMNS].to_csv(fh, sep="\t", index=False, header=False)
        summaries.append(chunk_out)
        print(f"  chunk {start}-{stop}: {n_scored_chunk} scored, {dropped_chunk} dropped "
              f"(<2 informative groups); t={time.time()-t0:.0f}s", flush=True)

    summary = pd.concat(summaries, ignore_index=True)
    # Validation required by the task: every scored variant has >= 2 informative groups.
    scored_mask = summary["z_score"].notna()
    assert (summary.loc[scored_mask, "n_groups"] >= 2).all(), "scored locus with <2 groups"
    print(f"n_scored={int(scored_mask.sum())}; "
          f"n_dropped_lt2_informative_groups={n_dropped}", flush=True)

    # Nulls for both statistics (matched controls).
    for score_col, suffix in (("replication_score", ""), ("subset_score", "_subset")):
        idx, p_emp, n_ctrl, bh = empirical_null(summary, score_col, args)
        summary.loc[idx, "n_controls" + suffix] = n_ctrl
        summary.loc[idx, "p_empirical" + suffix] = p_emp
        summary.loc[idx, "bh_fdr" + suffix] = bh
    summary = summary[COLUMNS]
    summary.to_csv(out_tsv, sep="\t", index=False)

    # ------------------------------------------------------------------ report #
    def quad_counts(p_col):
        d = summary.copy()
        d["absX"] = d["akbari_X"].abs()
        hi = d["absX"] >= SIG_ABS_X
        rep = d[p_col] <= 0.05
        return {
            "Q1_hi_hi": int((hi & rep).sum()),
            "Q2_hi_lo": int((hi & ~rep).sum()),
            "Q3_lo_hi": int((~hi & rep).sum()),
            "Q4_lo_lo": int((~hi & ~rep).sum()),
        }

    d = summary.copy()
    d["absX"] = d["akbari_X"].abs()
    finite_rep = np.isfinite(d["replication_score"])
    lo = d[finite_rep & (d["absX"] < SIG_ABS_X)]
    hi = d[finite_rep & (d["absX"] >= SIG_ABS_X)]
    # For the low-|X| pool, leave-one-out matched P (quadrant.empirical_p_leave_one_out).
    import quadrant
    d["p_use"] = d["p_empirical"]
    miss_lo = d["p_use"].isna() & (d["absX"] < SIG_ABS_X)
    if miss_lo.any():
        d.loc[miss_lo, "p_use_loo"] = quadrant.empirical_p_leave_one_out(d[miss_lo].copy())
        d["p_use"] = d["p_use"].fillna(d["p_use_loo"])
    d["rep_high"] = d["p_use"] <= 0.05
    lowfrac = float(d.loc[d["absX"] < SIG_ABS_X, "rep_high"].mean())
    highfrac = float(d.loc[d["absX"] >= SIG_ABS_X, "rep_high"].mean())

    fin = summary[np.isfinite(summary["subset_max_stat"])]
    pass_fdr = fin[fin["is_candidate"] & (fin["bh_fdr_subset"] <= 0.05)]
    prop_subset = pass_fdr["argmax_subset_size"] < pass_fdr["n_groups"]
    metrics = {
        "seed": SEED, "bin_years": BIN_YEARS,
        "sig_abs_x_threshold": SIG_ABS_X, "genome_wide_abs_x": GW_ABS_X,
        "n_qc_variants": int(pool.attrs["n_qc"]),
        "n_akbari_rows_scanned": int(pool.attrs["n_akbari_rows_scanned"]),
        "n_candidates_available": int(pool.attrs["n_cand_avail"]),
        "n_background_available": int(pool.attrs["n_bg_avail"]),
        "panel_cap_max_variants": args.max_variants,
        "n_panel_variants": int(n_panel),
        "n_panel_candidates": int(summary["is_candidate"].sum()),
        "n_panel_background": int((~summary["is_candidate"]).sum()),
        "n_trajectory_individuals": int(n_ind),
        "n_trajectory_groups": int(len(group_order)),
        "relatedness_matrix": REL_V3.name,
        "relatedness_units": "merged 6-group (trajectories_v3 merged_group)",
        "n_effective": float(len(group_order) ** 2 / r_full.sum()),
        "n_scored": int(scored_mask.sum()),
        "n_dropped_lt2_informative_groups": int(n_dropped),
        "null_method": "matched controls, same statistic (chrom, Akbari AF, "
                       "missingness, temporal coverage; same-chrom >=1Mb excluded; "
                       "calipers widen 1x/2x/4x)",
        "standard": {
            "quadrant_counts": quad_counts("p_empirical"),
            "frac_low_absX_p_le_0p05": lowfrac,
            "frac_high_absX_p_le_0p05": highfrac,
            "cand_bh_fdr_le_0p05": int((summary["is_candidate"] & (summary["bh_fdr"] <= 0.05)).sum()),
        },
        "max_over_subsets": {
            "n_candidates_finite_stat": int(len(fin[fin["is_candidate"]])),
            "n_candidates_bh_fdr_le_0p05": int(len(pass_fdr)),
            "n_pass_proper_subset": int(prop_subset.sum()),
            "candidate_subset_max_mean": float(fin.loc[fin["is_candidate"], "subset_max_stat"].mean()) if fin["is_candidate"].any() else None,
            "control_subset_max_mean": float(fin.loc[~fin["is_candidate"], "subset_max_stat"].mean()) if (~fin["is_candidate"]).any() else None,
        },
        "runtime_seconds": round(time.time() - t0, 1),
    }
    args.out_json.write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"wrote {out_tsv} ({len(summary)} rows) and {args.out_json}", flush=True)
    print(json.dumps({k: metrics[k] for k in (
        "n_panel_variants", "n_scored", "n_dropped_lt2_informative_groups",
        "runtime_seconds")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
