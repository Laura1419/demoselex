#!/usr/bin/env python3
"""Pre-registered subset-selection regime test (max-over-subsets statistic).

See ``docs/subset-regime-prereg.md`` (written before this script was run). This
is a targeted test of the ONE regime in which a within-trajectory / subset method
can beat Akbari's pooled GLMM: selection acting in a *subset* of population
histories while the rest stay flat or move the other way.

For each panel locus we reuse the per-trajectory weighted slopes and binomial SEs
from ``sim.trajectory_stat`` (SOURCES), then over every subset S of the
informative trajectory groups with 2 <= |S| <= 4 compute the correlated
fixed-effect GLS meta z restricted to S. The locus statistic is

    T = max over admissible S of |z(S)|,

calibrated against the SAME max-over-subsets statistic on matched control AADR
sites (caliper matcher of ``analysis.reanalysis_scan``). Candidate = Akbari
|X| >= 3.16. Decision rule: empirical upper-tail P then BH-FDR at alpha 0.05.

Writes ``data/derived/subset_regime.tsv`` incrementally after every genotype
chunk so a timeout loses nothing. Negative results are first-class.

Usage
-----
    .venv/bin/python analysis/subset_regime.py --max-sites 12000
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "analysis"))
sys.path.insert(0, str(REPO / "sim"))

import reanalysis_scan as rs  # noqa: E402
from trajectory_stat import _correlation, _fit_slopes  # noqa: E402

SIG_ABS_X = 3.16
BIN_YEARS = 1000
SEED = 42
MIN_SUBSET = 2
MAX_SUBSET = 4
N_TRAJ = 13  # trajectories 0..12 in data/derived/trajectories.tsv

# Precomputed admissible subset index tuples over the 13 groups (by size).
SUBSETS_BY_SIZE = {s: list(itertools.combinations(range(N_TRAJ), s))
                   for s in range(MIN_SUBSET, MAX_SUBSET + 1)}
ALL_SUBSETS = [c for s in range(MIN_SUBSET, MAX_SUBSET + 1)
               for c in SUBSETS_BY_SIZE[s]]

ESTIMATOR_COLUMNS = [
    "locus", "chrom", "pos", "akbari_X", "akbari_AF", "missingness",
    "temporal_coverage", "n_groups", "n_subsets", "subset_max_stat",
    "argmax_subset_size", "argmax_subset", "argmax_z", "complement_slope",
    "complement_z", "complement_flat", "complement_opposite", "is_candidate",
]
FINAL_COLUMNS = ESTIMATOR_COLUMNS + ["n_controls", "p_empirical", "bh_fdr"]


def gls_z(beta, se, r_sub):
    """Correlated fixed-effect GLS meta z for one subset (V=diag(se) R diag(se))."""
    v = se[:, None] * r_sub * se[None, :]
    w = np.linalg.inv(v)
    one = np.ones(len(beta))
    denom = float(one @ w @ one)
    meta = float(one @ w @ beta) / denom
    se_meta = float(np.sqrt(1.0 / denom))
    return meta, meta / se_meta


def subset_statistic(slopes_locus, r_full, traj_index):
    """Return T and argmax/complement diagnostics for one locus."""
    present = sorted(slopes_locus)
    present_set = set(present)
    idx = {t: traj_index[t] for t in present}
    slope = {t: slopes_locus[t][0] for t in present}
    se = {t: slopes_locus[t][1] for t in present}

    best = None
    n_sub = 0
    for combo in ALL_SUBSETS:
        if not set(combo).issubset(present_set):
            continue
        n_sub += 1
        labels = [f"t{c}" for c in combo]
        ridx = [idx[c] for c in combo]
        r_sub = r_full[np.ix_(ridx, ridx)]
        b = np.array([slope[c] for c in combo], float)
        s = np.array([se[c] for c in combo], float)
        meta, z = gls_z(b, s, r_sub)
        if best is None or abs(z) > abs(best[2]):
            best = (combo, meta, z)

    if best is None:
        return None
    combo, meta, z = best
    comp = tuple(t for t in present if t not in combo)
    comp_slope = comp_z = np.nan
    comp_flat = comp_opp = False
    if len(comp) >= 2:
        ridx = [idx[c] for c in comp]
        r_sub = r_full[np.ix_(ridx, ridx)]
        b = np.array([slope[c] for c in comp], float)
        s = np.array([se[c] for c in comp], float)
        c_meta, c_z = gls_z(b, s, r_sub)
        comp_slope, comp_z = c_meta, c_z
        comp_flat = abs(c_z) < 1.0
        comp_opp = np.sign(c_meta) != np.sign(meta) and c_meta != 0 and meta != 0
    return {
        "n_groups": len(present), "n_subsets": n_sub,
        "subset_max_stat": abs(z), "argmax_subset_size": len(combo),
        "argmax_subset": ",".join(f"t{c}" for c in combo), "argmax_z": z,
        "complement_slope": comp_slope, "complement_z": comp_z,
        "complement_flat": bool(comp_flat), "complement_opposite": bool(comp_opp),
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--max-sites", type=int, default=12000,
                    help="total cap on candidate+background sites analysed")
    ap.add_argument("--chunk-sites", type=int, default=2000)
    ap.add_argument("--min-controls", type=int, default=20)
    ap.add_argument("--af-caliper", type=float, default=0.10)
    ap.add_argument("--miss-caliper", type=float, default=0.10)
    ap.add_argument("--cov-caliper", type=float, default=0.15)
    ap.add_argument("--out-tsv", type=Path,
                    default=REPO / "data/derived/subset_regime.tsv")
    ap.add_argument("--out-json", type=Path,
                    default=REPO / "data/derived/subset_regime_metrics.json")
    args = ap.parse_args(argv)

    t0 = time.time()
    rel_path = REPO / "data/derived/traj_relatedness.tsv"
    rel = None
    relatedness_used = rel_path.exists()
    if relatedness_used:
        rel = pd.read_csv(rel_path, sep="\t", index_col=0)
        rel.index = [f"t{i}" for i in rel.index]
        rel.columns = [f"t{i}" for i in rel.columns]
    traj_index = {i: i for i in range(N_TRAJ)}
    r_full = (_correlation(rel, [f"t{i}" for i in range(N_TRAJ)])
              if rel is not None else np.eye(N_TRAJ))
    print(f"relatedness applied={relatedness_used}; groups={N_TRAJ}; "
          f"subsets per locus <= {len(ALL_SUBSETS)}", flush=True)

    cohort, n_unmapped = rs.load_cohort(REPO / "data/derived/trajectories.tsv")
    rows = np.unique(cohort["row"].to_numpy())
    row_pos = {int(r): k for k, r in enumerate(rows)}
    cohort_pos = cohort["row"].map(row_pos).to_numpy(np.int64)
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

    print(f"cohort: {len(cohort)} individuals, {n_groups} (trajectory,bin) cells; "
          f"unmapped={n_unmapped}", flush=True)

    panel_info = rs.build_panel(args.max_sites)
    panel = panel_info["panel"]
    is_cand = panel_info["is_cand"]
    labels = panel_info["labels"]
    n_panel = len(panel)
    n_ind = len(rows)

    gm = rs.GenotypeMatrix(str(REPO / rs.BUNDLE / (rs.MAIN_PREFIX + ".geno")))
    assert gm.n_ind == rs.parse_header(gm.path).n_ind

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
        codes = rs.read_rows_sites(gm, rows, sites)
        alt_mat = np.where(codes <= 2, codes, 0).astype(np.int32)
        tot_mat = (2 * (codes <= 2)).astype(np.int32)
        alt_sum = np.add.reduceat(alt_mat[sorter], starts, axis=0)
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
        pair = tidy["locus"] + "\x1f" + tidy["trajectory"]
        per = tidy.assign(pair=pair).groupby("pair", sort=False)["time"].nunique()
        keep_pairs = set(per[per >= 2].index)
        tidy = tidy[pair.isin(keep_pairs)]
        ntraj = tidy.groupby("locus", sort=False)["trajectory"].nunique()
        good_loci = set(ntraj[ntraj >= 2].index)

        slopes = _fit_slopes(tidy)
        # Per-locus dict: trajectory index -> (slope, se).
        sl_by_locus: dict[str, dict[int, tuple[float, float]]] = {}
        for row in slopes.itertuples(index=False):
            ti = int(row.trajectory[1:])
            sl_by_locus.setdefault(row.locus, {})[ti] = (row.slope, row.se)

        chunk_rows = []
        for k in range(stop - start):
            lab = labels[start + k]
            rec = {"locus": lab, "chrom": panel_info["chrom"][start + k],
                   "pos": panel_info["pos"][start + k],
                   "akbari_X": panel_info["akbari_X"][start + k],
                   "akbari_AF": panel_info["akbari_AF"][start + k],
                   "missingness": miss[start + k],
                   "temporal_coverage": tcov[start + k],
                   "is_candidate": bool(is_cand[start + k])}
            sl = sl_by_locus.get(lab)
            if sl is None or len(sl) < MIN_SUBSET:
                rec.update({c: np.nan for c in ESTIMATOR_COLUMNS[7:17]})
                rec["n_groups"] = len(sl) if sl else 0
                rec["n_subsets"] = 0
                chunk_rows.append(rec)
                continue
            st = subset_statistic(sl, r_full, traj_index)
            rec.update(st)
            chunk_rows.append(rec)
        n_dropped_support += (stop - start) - len(good_loci)
        chunk_df = pd.DataFrame(chunk_rows)
        with open(out_tsv, "a", encoding="utf-8") as fh:
            chunk_df[ESTIMATOR_COLUMNS].to_csv(fh, sep="\t", index=False,
                                               header=False)
        summaries.append(chunk_df)
        print(f"  chunk {start}-{stop}: {len(chunk_df)} loci, "
              f"dropped={ (stop-start)-len(good_loci)}; t={time.time()-t0:.0f}s",
              flush=True)

    summary = pd.concat(summaries, ignore_index=True)
    finite = np.isfinite(summary["subset_max_stat"].to_numpy(dtype=float))
    for col in ("n_controls", "p_empirical", "bh_fdr"):
        summary[col] = np.nan
    cand = summary[summary["is_candidate"] & finite].copy()
    bg = summary[(~summary["is_candidate"]) & finite].copy()

    p_emp = np.full(len(cand), np.nan)
    n_ctrl = np.zeros(len(cand), dtype=int)
    for j, (_cidx, crow) in enumerate(cand.iterrows()):
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
                break
        if chosen is None:
            chosen = pool
        m = len(chosen)
        n_ctrl[j] = m
        if m == 0:
            continue
        scores = chosen["subset_max_stat"].to_numpy(dtype=float)
        cs = float(crow["subset_max_stat"])
        p_emp[j] = (1 + int((scores >= cs).sum())) / (m + 1)

    cand["p_empirical"] = p_emp
    cand["n_controls"] = n_ctrl
    cand["bh_fdr"] = np.nan
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
        cand.loc[ok, "bh_fdr"] = out

    summary.loc[cand.index, "n_controls"] = n_ctrl
    summary.loc[cand.index, "p_empirical"] = p_emp
    summary.loc[cand.index, "bh_fdr"] = cand["bh_fdr"].to_numpy()
    summary = summary[FINAL_COLUMNS]
    summary.to_csv(out_tsv, sep="\t", index=False)

    cand_fin = summary.loc[summary["is_candidate"] & finite]
    pass_fdr = cand_fin[cand_fin["bh_fdr"] <= 0.05]
    lowx = pass_fdr[pass_fdr["akbari_X"].abs() < SIG_ABS_X]
    prop = pass_fdr["argmax_subset_size"] < pass_fdr["n_groups"]
    cand_mean = float(cand_fin["subset_max_stat"].mean()) if len(cand_fin) else None
    ctrl_mean = float(bg["subset_max_stat"].mean()) if len(bg) else None
    metrics = dict(panel_info["stats"], **{
        "seed": SEED, "min_subset": MIN_SUBSET, "max_subset": MAX_SUBSET,
        "alpha": 0.05, "sig_abs_x_threshold": SIG_ABS_X,
        "n_loci_scored": int(len(summary)),
        "n_loci_finite_stat": int(finite.sum()),
        "n_candidates_scored": int(summary["is_candidate"].sum()),
        "n_background_scored": int((~summary["is_candidate"]).sum()),
        "n_candidates_finite_p": int(np.isfinite(p_emp).sum()),
        "n_candidates_bh_fdr_le_0p05": int((cand_fin["bh_fdr"] <= 0.05).sum()),
        "n_pass_low_absX": int(len(lowx)),
        "n_pass_proper_subset": int(prop.sum()),
        "candidate_subset_max_mean": cand_mean,
        "control_subset_max_mean": ctrl_mean,
        "relatedness_applied": bool(relatedness_used),
        "null_method": "matched controls, same max-over-subsets statistic "
                       "(chrom, AF, missingness, temporal coverage; same-chrom "
                       ">=1Mb excluded; calipers widen 1x/2x/4x)",
        "runtime_seconds": round(time.time() - t0, 1),
    })
    args.out_json.write_text(json.dumps(metrics, indent=2) + "\n")
    print(f"wrote {out_tsv} ({len(summary)} rows) and {args.out_json}", flush=True)
    print(json.dumps({k: metrics[k] for k in (
        "n_loci_scored", "n_candidates_scored", "n_candidates_finite_p",
        "n_candidates_bh_fdr_le_0p05", "n_pass_low_absX",
        "candidate_subset_max_mean", "control_subset_max_mean",
        "runtime_seconds")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
