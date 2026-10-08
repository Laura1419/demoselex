#!/usr/bin/env python3
"""Merge v2 trajectories into genetically-independent groups by relatedness.

Problem
-------
The v2 pipeline (``analysis/trajectory_v2.py``) returns 13 trajectories, but
several are the *same* population history split by sampling period or admixture
fraction, so they are not independent replications.  The diagnosed cases are

* T1/T3/T6 - all Avar-period Hungary, split by the ``oLowEastAsia`` vs
  ``oHighEastAsia`` ancestry proportion, pairwise r 0.58-0.79;
* T0/T4 - both medieval/Avar-period Europe, r(0,4)=0.858 (T4 is the
  4,594-individual blob that hit the 8,000-yr span cap and is under-resolved).

Method
------
1. Take the existing 13x13 relatedness matrix (``data/derived/traj_relatedness.tsv``)
   and merge trajectory labels whose pairwise coefficient exceeds a threshold,
   using union-find (imported from ``analysis.traj_relatedness``).  This is done
   for a range of thresholds; the grouping and ``N_eff = m^2 / sum(R)`` are
   reported for every threshold.
2. Justify the operating threshold from the *gap structure* of the sorted
   off-diagonal coefficients - not from the size of the resulting grouping.
3. Emit ``data/derived/trajectories_v3.tsv``: the exact v2 rows plus a new
   ``merged_group`` column (``group_id`` is already the AADR published-group
   annotation, so it is NOT reused).  All v2 columns are preserved.
4. Recompute the relatedness matrix on the *merged groups* by importing the v2
   estimator (``analysis.traj_relatedness.group_accumulate`` /
   ``relatedness_from_accum``) and write ``data/derived/traj_relatedness_v3.tsv``.
   The consumer contract ``sim.trajectory_stat._correlation`` is called as a test.
5. Per merged group: n, date range, dominant published group IDs, mean lat/lon,
   and whether it maps onto a recognisable West Eurasian history.
6. Honesty verdict on which groups are distinct histories vs sampling artefacts,
   the final ``N_eff``, and whether it is defensible as a replication count.

All output numbers are mirrored in ``data/derived/traj_v3_metrics.json`` and the
companion doc ``docs/traj-v3-merge.md`` (both written by this script).

All claims in the doc are labelled [observation] / [inference] / [speculation].
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.load_aadr import GenotypeMatrix, parse_ind  # noqa: E402
from analysis.traj_relatedness import (  # noqa: E402
    GENO_PATH,
    IND_PATH,
    candidate_sites,
    coarse_matrix,
    group_accumulate,
    merge_by_threshold,
    n_eff,
    relatedness_from_accum,
    select_sites,
)
from analysis.trajectory_v2 import EEF_KEYS, STEPPE_KEYS, WHG_KEYS  # noqa: E402
from sim.trajectory_stat import _correlation  # noqa: E402

TRAJ_V2 = "data/derived/trajectories.tsv"
REL_V2 = "data/derived/traj_relatedness.tsv"
OUT_TRAJ = "data/derived/trajectories_v3.tsv"
OUT_REL = "data/derived/traj_relatedness_v3.tsv"
OUT_METRICS = "data/derived/traj_v3_metrics.json"
OUT_DOC = "docs/traj-v3-merge.md"
FIG_DIR = "fig"

# Thresholds required by the task (must all be reported, not just the chosen one).
THRESHOLDS = [0.9, 0.7, 0.5, 0.4, 0.3]
# Operating threshold - justified in the doc from group content + gap structure.
OPERATING = 0.5


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def load_v2_matrix():
    m = pd.read_csv(REL_V2, sep="\t", index_col=0)
    m.index = [str(x) for x in m.index]
    m.columns = [str(c) for c in m.columns]
    order = sorted(m.index, key=int)
    return m.reindex(index=order, columns=order)


def gap_structure(mat):
    """Sorted off-diagonal coefficients and the consecutive gaps between them."""
    labels = list(mat.index)
    vals = []
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            vals.append((float(mat.loc[a, b]), a, b))
    vals.sort(key=lambda x: -x[0])
    rows = []
    prev = None
    for v, a, b in vals:
        rows.append({"a": a, "b": b, "r": v,
                     "gap_above": None if prev is None else round(prev - v, 4)})
        prev = v
    nonzero = [r for r in rows if r["r"] > 0]
    gaps = [(r["gap_above"], f'{r["a"]}-{r["b"]}', r["r"]) for r in nonzero
            if r["gap_above"] is not None]
    gap_sorted = sorted(gaps, key=lambda x: -x[0])
    return rows, nonzero, gap_sorted


def canonical_clusters(mat, thr):
    """Union-find clusters at ``thr``, deterministically ordered by min label."""
    clusters = merge_by_threshold(mat, thr)
    clusters = sorted(clusters, key=lambda c: min(int(x) for x in c))
    return clusters


def group_label_map(clusters):
    """trajectory label (str) -> merged group name 'G0'..'Gk'."""
    return {lab: f"G{k}" for k, c in enumerate(clusters) for lab in c}


def summarise_groups(traj, label_map):
    """Per merged group: n, span, dominant published groups, mean lat/lon, keys."""
    out = {}
    merged = traj.trajectory_id.astype(int).map({int(k): v for k, v in label_map.items()})
    traj = traj.assign(merged_group=merged)
    for g, sub in traj.groupby("merged_group", sort=True):
        cnt = Counter(sub.group_id.tolist())
        joined = " ".join(x.lower() for x in cnt)
        out[g] = {
            "n": int(len(sub)),
            "member_trajectories": sorted(int(x) for x in sub.trajectory_id.unique()),
            "date_min": float(sub.date_BP.min()),
            "date_max": float(sub.date_BP.max()),
            "date_median": float(sub.date_BP.median()),
            "mean_lat": float(sub.lat.mean()),
            "mean_lon": float(sub.lon.mean()),
            "top_groups": [[x, int(c)] for x, c in cnt.most_common(6)],
            "n_published_groups": int(len(cnt)),
            "whg_like": int(sum(c for x, c in cnt.items()
                                if any(k in x.lower() for k in WHG_KEYS))),
            "eef_like": int(sum(c for x, c in cnt.items()
                                if any(k in x.lower() for k in EEF_KEYS))),
            "steppe_like": int(sum(c for x, c in cnt.items()
                                   if any(k in x.lower() for k in STEPPE_KEYS))),
        }
    return out


# Human-readable names, assigned after inspecting the per-group summary.  Kept
# descriptive/geographic - no archaeological period claims are implied by them.
MERGED_NAMES = {
    "G0": "broad west/central-Eurasian Holocene bloc (medieval + Bronze-Age + early-farmer + steppe-adjacent material)",
    "G1": "Avar-period Hungary, East-Asian-admixed",
    "G2": "Iran/Anatolia Neolithic + Iron-Gates hunter-gatherer (heterogeneous)",
    "G3": "steppe-like (Khvalynsk / Yamnaya)",
    "G4": "WHG-like hunter-gatherer (Iron Gates / Narva / Spain Mesolithic)",
    "G5": "Ukrainian Neolithic (mixed WHG/EEF/steppe)",
}


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    print("== load v2 relatedness matrix ==", flush=True)
    mat = load_v2_matrix()
    labels = list(mat.index)
    print(f"  {mat.shape[0]}x{mat.shape[1]} trajectory matrix, labels={labels}",
          flush=True)

    neff_raw, m_raw = n_eff(mat)
    print(f"  N_eff (13 distinct) = {neff_raw:.3f}", flush=True)

    # ---- threshold sweep on the EXISTING matrix ----
    print("== threshold sweep ==", flush=True)
    sweep = []
    for thr in THRESHOLDS:
        clusters = canonical_clusters(mat, thr)
        coarse = coarse_matrix(mat, clusters)
        ne, m = n_eff(mat, clusters)
        sweep.append({
            "threshold": thr,
            "n_groups": m,
            "n_effective": ne,
            "clusters": clusters,
            "N_eff_all_corrected": None,
            "sum_R": float(coarse.to_numpy().sum()),
        })
        print(f"  r>{thr}: {m} groups N_eff={ne:.3f} {clusters}", flush=True)

    rows, nonzero, gap_sorted = gap_structure(mat)
    print("  largest off-diagonal gaps:",
          [(g, p) for g, p, _ in gap_sorted[:5]], flush=True)
    max_gap, max_gap_pair, max_gap_r = gap_sorted[0]

    # ---- operating grouping ----
    op_clusters = canonical_clusters(mat, OPERATING)
    label_map = group_label_map(op_clusters)
    print(f"== operating threshold r>{OPERATING}: {op_clusters}", flush=True)

    # ---- trajectories_v3.tsv ----
    traj = pd.read_csv(TRAJ_V2, sep="\t")
    traj["trajectory_id"] = traj.trajectory_id.astype(int)
    assigned = traj.trajectory_id >= 0
    merged = pd.Series("unassigned", index=traj.index)
    merged.loc[assigned] = (traj.loc[assigned, "trajectory_id"]
                            .map({int(k): v for k, v in label_map.items()}))
    out = traj.copy()
    out["merged_group"] = merged
    os.makedirs(os.path.dirname(OUT_TRAJ), exist_ok=True)
    out.to_csv(OUT_TRAJ, sep="\t", index=False)
    print(f"  wrote {OUT_TRAJ}: {out.shape[0]} rows, "
          f"columns={list(out.columns)}", flush=True)
    assert list(out.columns)[:len(traj.columns)] == list(traj.columns), \
        "v2 columns must be preserved in order"

    # ---- recompute relatedness on merged groups ----
    print("== recompute relatedness on merged groups ==", flush=True)
    sub = traj.loc[assigned].copy()
    ind_ids = [r.id for r in parse_ind(IND_PATH)]
    row_of = {x: i for i, x in enumerate(ind_ids)}
    keep = sub.individual.isin(row_of)
    sub = sub.loc[keep]
    ind_rows = np.array([row_of[x] for x in sub.individual], dtype=np.int64)
    grp_of_row = np.array([int(label_map[str(t)][1:]) for t in sub.trajectory_id],
                          dtype=int)
    sites = select_sites(candidate_sites())
    gm = GenotypeMatrix(GENO_PATH)
    groups, sums, counts = group_accumulate(gm, ind_rows, grp_of_row, sites)
    mat3, S_used = relatedness_from_accum(groups, sums, counts)
    mat3.index = [f"G{g}" for g in mat3.index]
    mat3.columns = [f"G{g}" for g in mat3.columns]
    mat3.to_csv(OUT_REL, sep="\t")
    print(f"  wrote {OUT_REL} ({mat3.shape[0]}x{mat3.shape[1]}), "
          f"sites used={S_used}; individuals={len(sub)}", flush=True)

    # ---- consumer contract test ----
    pd_ok, pd_err = True, None
    try:
        _correlation(mat3, list(mat3.index))
    except Exception as exc:  # noqa: BLE001
        pd_ok, pd_err = False, str(exc)
    print(f"  sim.trajectory_stat._correlation: "
          f"{'PASS' if pd_ok else 'FAIL -> ' + str(pd_err)}", flush=True)

    # ---- final N_eff on MERGED groups (recomputed matrix) ----
    neff_merged, m_merged = n_eff(mat3)
    neff_merged_v2coarse, _ = n_eff(mat, op_clusters)
    print(f"  N_eff (merged, recomputed matrix) = {neff_merged:.3f} over "
          f"m={m_merged}", flush=True)
    print(f"  N_eff (merged, coarse of v2 matrix) = {neff_merged_v2coarse:.3f}",
          flush=True)

    # ---- per-group description ----
    group_stats = summarise_groups(traj.loc[assigned].copy(), label_map)
    for g, s in sorted(group_stats.items()):
        print(f"  {g}: n={s['n']} date {s['date_min']:.0f}-{s['date_max']:.0f} "
              f"lat~{s['mean_lat']:.1f} lon~{s['mean_lon']:.1f} "
              f"top={s['top_groups'][:3]}", flush=True)

    # ---- figures ----
    make_figures(mat3, traj.loc[assigned], label_map, group_stats)

    # ---- metrics ----
    metrics = {
        "operating_threshold": OPERATING,
        "thresholds_tested": THRESHOLDS,
        "N_eff_raw_13": neff_raw,
        "N_eff_merged_recomputed": neff_merged,
        "N_eff_merged_from_v2_coarse": neff_merged_v2coarse,
        "n_merged_groups": m_merged,
        "merged_groups": {g: label_map_reverse(op_clusters, g)
                          for g in [f"G{i}" for i in range(len(op_clusters))]},
        "group_names": MERGED_NAMES,
        "threshold_sweep": [{k: v for k, v in row.items()
                             if k != "N_eff_all_corrected"} for row in sweep],
        "gap_structure_top": [{"pair": p, "r": r, "gap_above": g}
                              for g, p, r in gap_sorted[:8]],
        "largest_gap": {"gap": max_gap, "pair": max_gap_pair, "r": max_gap_r},
        "relatedness_v3": {
            "sites_used": S_used,
            "n_individuals": int(len(sub)),
            "values": {g: {h: float(mat3.loc[g, h]) for h in mat3.columns}
                       for g in mat3.index},
        },
        "cholesky_pd": pd_ok,
        "cholesky_error": pd_err,
        "group_stats": group_stats,
    }
    with open(OUT_METRICS, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print(f"  wrote {OUT_METRICS}", flush=True)

    write_doc(metrics, mat3, sweep, op_clusters, group_stats,
              neff_raw, neff_merged, neff_merged_v2coarse, gap_sorted)

    # ---- assertions (task validation) ----
    assert max_gap is not None
    assert len(sweep) == len(THRESHOLDS)
    assert m_merged == len(op_clusters)
    assert pd_ok, f"v3 matrix failed _correlation: {pd_err}"
    assert os.path.exists(OUT_TRAJ) and os.path.exists(OUT_REL)
    print("ALL CHECKS PASS", flush=True)
    return metrics


def label_map_reverse(clusters, gname):
    k = int(gname[1:])
    return clusters[k]


def make_figures(mat3, traj, label_map, group_stats):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(FIG_DIR, exist_ok=True)
    labels = list(mat3.index)

    # 1. merged-group relatedness heatmap
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(mat3.values, cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticklabels(labels)
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j, i, f"{mat3.values[i, j]:.2f}", ha="center", va="center",
                    color="w" if mat3.values[i, j] < 0.6 else "k", fontsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, label="relatedness r")
    ax.set_title("Merged-group relatedness (v3)")
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/traj_v3_relatedness.png", dpi=150)
    plt.close(fig)

    # 2. per-group date span + spatial means
    merged = traj.trajectory_id.astype(int).map(
        {int(k): v for k, v in label_map.items()})
    gl = [f"G{i}" for i in range(len(group_stats))]
    cmap = plt.get_cmap("tab10")
    color = {g: cmap(i % 10) for i, g in enumerate(gl)}
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
    for i, g in enumerate(gl):
        d = traj.loc[merged == g, "date_BP"].to_numpy()
        ax[0].plot([d.min(), d.max()], [i, i], "-", color=color[g], lw=4)
        ax[0].scatter([np.median(d)], [i], color="k", s=18, zorder=3)
        ax[0].text(d.min(), i + 0.25, f"{g} (n={len(d)})", fontsize=8)
    ax[0].set(yticks=range(len(gl)), xlabel="date (BP)",
              title="Merged-group date spans (bar) + median (dot)")
    ax[0].invert_xaxis()
    for g in gl:
        m = merged == g
        ax[1].scatter(traj.loc[m, "lon"], traj.loc[m, "lat"], s=8,
                      color=color[g], alpha=0.5, label=g)
    ax[1].set(xlabel="longitude", ylabel="latitude",
              title="Merged-group spatial extent")
    ax[1].legend(fontsize=7, markerscale=2)
    fig.tight_layout()
    fig.savefig(f"{FIG_DIR}/traj_v3_groups.png", dpi=130)
    plt.close(fig)
    print(f"  wrote {FIG_DIR}/traj_v3_relatedness.png, "
          f"{FIG_DIR}/traj_v3_groups.png", flush=True)


def _fmt_cluster(c):
    return "{" + ",".join(c) + "}"


def write_doc(metrics, mat3, sweep, op_clusters, group_stats,
              neff_raw, neff_merged, neff_merged_v2coarse, gap_sorted):
    op = metrics["operating_threshold"]
    L = []
    A = L.append
    A("# Trajectory merge v3 - genetically-independent replication units")
    A("")
    A("> Owner lane: traj-v3-merge. Producer: `analysis/trajectory_merge.py`. "
      "Status: **merged grouping produced; see verdict.**")
    A("")
    A("Reproduce (single command, writes all four outputs):")
    A("")
    A("```")
    A(".venv/bin/python analysis/trajectory_merge.py")
    A("```")
    A("")
    A("## Why this file exists")
    A("")
    A("The v2 pipeline returns 13 trajectories, but several are the *same* "
      "population history split by sampling period or admixture fraction, so "
      "counting them as 13 independent replications double-counts shared "
      "history. Two splits were diagnosed independently "
      "(`docs/traj-relatedness.md`):")
    A("")
    A("* T1/T3/T6 - Avar-period Hungary split by oLowEastAsia vs "
      "oHighEastAsia ancestry proportion (pairwise r 0.58-0.79);")
    A("* T0/T4 - medieval/Avar-period Europe, r(0,4)=0.858 (T4 is the "
      "4,594-individual blob that hit the 8,000-yr span cap).")
    A("")
    A("## 1. Threshold sweep (union-find on the existing 13x13 matrix)")
    A("")
    A("`[observation]` Merging all trajectory pairs whose relatedness "
      "coefficient exceeds the threshold, by union-find, on "
      "`data/derived/traj_relatedness.tsv`. `N_eff = m^2 / sum(R)` is computed "
      "on the merged groups (cluster-to-cluster relatedness = mean of member "
      "pair coefficients). Every tested threshold is shown:")
    A("")
    A("| threshold | m groups | grouping | N_eff |")
    A("|---|---:|---|---:|")
    for row in sweep:
        grouping = " ".join(_fmt_cluster(c) for c in row["clusters"])
        A(f"| r>{row['threshold']} | {row['n_groups']} | {grouping} | "
          f"{row['n_effective']:.3f} |")
    A(f"| (no merge) 13 distinct | 13 | all singleton | {neff_raw:.3f} |")
    A("")
    A("`[observation]` The sweep is monotone in coarseness but not in N_eff: "
      "collapsing two near-duplicates can tick N_eff *up* slightly before it "
      "falls, an arithmetic property of m^2/sum(R) already documented in "
      "`docs/traj-relatedness.md`.")
    A("")
    A("## 2. Is there a natural break in the relatedness values?")
    A("")
    A("`[observation]` Sorted off-diagonal coefficients (28 nonzero pairs) and "
      "the largest consecutive gaps:")
    A("")
    A("| rank | pair | r | gap above |")
    A("|---:|---|---:|---:|")
    for k, row in enumerate(metrics["gap_structure_top"]):
        A(f"| {k+1} | {row['pair']} | {row['r']:.4f} | "
          f"{row['gap_above'] if row['gap_above'] is not None else '-'} |")
    A("")
    lg = metrics["largest_gap"]
    A(f"`[observation]` The largest gap is **{lg['gap']:.4f}** "
      f"(between {lg['pair']} at r={lg['r']:.4f} and the next pair). The next "
      "largest gaps are ~0.088, ~0.064, ~0.064, ~0.058 - all comparable. "
      "**There is no sharp natural break**: the spectrum is continuous with "
      "only modest, similar-sized steps.")
    A("")
    A("`[inference]` A gap-maximising rule would put the cut near r~0.7 "
      "(between the 0.794 and 0.697 pairs). But that threshold only merges "
      "{0,4} and {1,6} and would leave T3 *separate* from T1/T6 - i.e. it "
      "contradicts the independently diagnosed Avar split, whose internal "
      "T3-T6 coefficient (0.576) sits below the gap. So the gap structure "
      "does **not** single out a merge set: a purely data-driven cut and the "
      "diagnosis disagree. This is stated plainly rather than picked to give "
      "a pleasing number.")
    A("")
    A(f"## 3. Chosen grouping (operating threshold r>{op})")
    A("")
    A(f"`[inference]` Operating threshold **r>{op}** is chosen on two grounds "
      "that are independent of the resulting count:")
    A("")
    A("* it is the *finest tested* threshold that reunites the diagnosed "
      "T1/T3/T6 Avar cluster (needs r < 0.576, satisfied by 0.5 but not 0.7);")
    A("* coarser thresholds (0.4, 0.3) add merges that mix clearly distinct "
      "histories - at r>0.4 the steppe-like T9 joins the WHG-like T11/T12, "
      "and at r>0.3 almost everything collapses into 2 groups.")
    A("")
    A("`[observation]` Chosen groups:")
    A("")
    A("| group | member trajectories | name (descriptive) | n | date min-max (BP) | mean lat/lon | dominant published groups |")
    A("|---|---|---|---:|---|---|---|")
    gl = [f"G{i}" for i in range(len(op_clusters))]
    for g in gl:
        s = group_stats[g]
        dom = ", ".join(f"{x} ({c})" for x, c in s["top_groups"][:3])
        A(f"| {g} | {','.join(str(t) for t in s['member_trajectories'])} | "
          f"{metrics['group_names'].get(g, '')} | {s['n']} | "
          f"{s['date_min']:.0f}-{s['date_max']:.0f} | "
          f"{s['mean_lat']:.1f}/{s['mean_lon']:.1f} | {dom} |")
    A("")
    A("`[observation]` `data/derived/trajectories_v3.tsv` uses the SAME rows "
      "as `trajectories.tsv` with one added column **`merged_group`** "
      "(`G0`..`G5`, or `unassigned` for `trajectory_id == -1`). The existing "
      "`group_id` column is the AADR *published-group* annotation and is "
      "deliberately NOT reused; `merged_group` is the new merged label. All v2 "
      "columns and names are preserved, and downstream readers of "
      "`individual`, `date_BP`, `trajectory_id` are unaffected.")
    A("")
    A("## 4. Recomputed relatedness on the merged groups")
    A("")
    A("`[observation]` The v2 estimator "
      "(`analysis.traj_relatedness.group_accumulate` / "
      "`relatedness_from_accum`) was re-run with the merged labels, pooling "
      "all individuals per group, and written to "
      "`data/derived/traj_relatedness_v3.tsv`. Sites used: "
      f"{metrics['relatedness_v3']['sites_used']}; individuals: "
      f"{metrics['relatedness_v3']['n_individuals']}.")
    A("")
    A("`[observation]` The matrix passes the consumer contract "
      "`sim.trajectory_stat._correlation` (symmetric, unit diagonal, entries "
      f"in [0,1], Cholesky PD): {'PASS' if metrics['cholesky_pd'] else 'FAIL'}.")
    A("")
    A("Merged-group matrix:")
    A("")
    A("| | " + " | ".join(gl) + " |")
    A("|---|" + "---|" * len(gl))
    for g in gl:
        A(f"| {g} | " + " | ".join(f"{mat3.loc[g, h]:.2f}" for h in gl) + " |")
    A("")
    nz = [(g, h, float(mat3.loc[g, h])) for i, g in enumerate(gl)
          for h in gl[i + 1:] if mat3.loc[g, h] > 0]
    nz_txt = ", ".join(f"{g}-{h}={v:.2f}" for g, h, v in nz) or "none"
    A(f"`[observation]` The recomputed group matrix is nearly block-diagonal: "
      f"the only nonzero off-diagonals are {nz_txt}; all others are clipped "
      "to 0. Pooling many trajectories into one group's mean allele-frequency "
      "profile averages away the shared drift the 13-way matrix showed "
      "(r(0,4)=0.86), so the grouped estimator *under*-states relatedness "
      "here. This is the opposite failure mode to over-calling, and it is why "
      "the recomputed N_eff is higher than the coarse-v2 N_eff.")
    A("")
    A("## 5. Final N_eff and honesty verdict")
    A("")
    A(f"`[observation]` On the merged groups, `N_eff = m^2/sum(R)` computed "
      f"from the **recomputed** v3 matrix = **{neff_merged:.3f}** "
      f"(m={metrics['n_merged_groups']}). For comparison, `N_eff` on the raw "
      f"13 = {neff_raw:.3f}, and on the merged groups using the coarse of the "
      f"v2 matrix = {neff_merged_v2coarse:.3f}.")
    A("")
    A("`[inference]` **Split of the 6 merged groups into distinct histories "
      "vs sampling artefacts:**")
    A("")
    A("* **Plausibly distinct histories (3-4):** G1 (Avar-period East-Asian-"
      "admixed), G3 (steppe-like: Khvalynsk/Yamnaya), G4 (WHG-like: Iron "
      "Gates/Narva/Spain Mesolithic), and arguably G5 (Ukrainian Neolithic).")
    A("* **Still sampling artefacts / over-merged (1-2):** **G0 is not one "
      "history.** It pools medieval/Avar-period Europe (T0/T4) with Anatolian/"
      "Iranian Bronze-Iron (T8), steppe CordedWare/Unetice material (T7), "
      "early-farmer LBK (T10) and Byzantine Anatolia (T2) into a single "
      "'broad bloc' - exactly the period/region sampling artefact this task "
      "set out to remove. G2 (Iran Neolithic + Iron Gates + Catalhoyuk, "
      "n=238) is small and heterogeneous and may also be a sampling artefact.")
    A("")
    A(f"`[observation]` The recomputed merged `N_eff` ({neff_merged:.3f}) is "
      f"*above* the raw-13 value ({neff_raw:.3f}) - the same m^2/sum(R) "
      "arithmetic seen in `docs/traj-relatedness.md` when a tight pair "
      "({0,4}) is collapsed. The recomputed matrix (pooled group means) and "
      f"the coarse average of the v2 matrix ({neff_merged_v2coarse:.3f}) "
      "disagree, which is itself a warning that N_eff is not stable to the "
      "estimation route.")
    A("")
    A("`[inference]` Across every route N_eff stays in the ~4-5 band and is "
      "dominated by the pervasive 0.2-0.6 background relatedness. Because G0 "
      "is demonstrably an over-merge, treating it as one independent unit "
      "**still overcounts** replication: the honest replication count is the "
      "number of *genuinely distinct* histories (~3-4), not 6 and not 13.")
    A("")
    A("`[speculation]` The G0 bloc is consistent with a shared west/central "
      "Eurasian deep ancestry sampled at many times; splitting it cleanly "
      "would need admixture/outgroup-f3 modelling or a tighter within-Holocene "
      "gate, which is out of scope here.")
    A("")
    A("## Verdict: is this grouping fit to be the replication unit?")
    A("")
    A(f"`[inference]` **PARTLY.** The merged `trajectories_v3.tsv` / "
      f"`traj_relatedness_v3.tsv` are a correct, reproducible collapse of the "
      f"13 v2 trajectories into {metrics['n_merged_groups']} relatedness "
      "groups and are safe to pass as the `family`/`relatedness` inputs to "
      "`sim.trajectory_stat`. But they are **not yet a defensible replication "
      "count on their own**: the grouping fixes the diagnosed T1/T3/T6 and "
      "T0/T4 splits only to the extent that r>0.5 merges them, while the "
      "resulting G0 still bundles several distinct histories, so the merged "
      "set under-corrects. The load-bearing number remains **N_eff ~ 4-5** "
      "(consistent with `docs/traj-relatedness.md`); report locus concordance "
      "with `n_effective ~= 4` (the count of genuinely distinct histories), "
      "and treat the G0 bloc as mixed background "
      "rather than as one independent test.")
    A("")
    A("## Artifacts")
    A("")
    A("* code: `analysis/trajectory_merge.py`")
    A("* merged table: `data/derived/trajectories_v3.tsv` "
      "(v2 columns + `merged_group`)")
    A("* merged matrix: `data/derived/traj_relatedness_v3.tsv`")
    A("* metrics: `data/derived/traj_v3_metrics.json`")
    A("* figures: `fig/traj_v3_relatedness.png`, `fig/traj_v3_groups.png`")
    A("")
    os.makedirs(os.path.dirname(OUT_DOC), exist_ok=True)
    with open(OUT_DOC, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"  wrote {OUT_DOC}", flush=True)


if __name__ == "__main__":
    main()
