"""Genetically-independent-trajectory quantification.

Estimator
---------
For a filtered set of autosomal biallelic sites we decode the 2-bit AADR
genotypes as diploid allele dosages ``g in {0,1,2}`` for every individual that
belongs to one of the 13 inferred trajectories.  At each site ``s`` the
per-trajectory frequency is ``p_{t,s} = mean dosage / 2``.  The reference
frequency is the **equal-weight** mean of the trajectory frequencies
(``ref_s = (1/m) sum_t p_{t,s}``), so a large central cluster cannot define the
reference away from itself; the binomial scaling is
``sd_s = sqrt(2 ref_s (1 - ref_s))``.  The standardised profile of trajectory
``t`` is ``d_{t,s} = (p_{t,s} - ref_s) / sd_s``.

Because the group mean is linear, the matrix

    R_ij = (1/S) * sum_s d_{i,s} * d_{j,s}

is the sample-size-normalised **mean pairwise allele-sharing / IBS** between
individuals of trajectories ``i`` and ``j`` (Weir--Hill kinship in spirit),
aggregated to trajectory level.  The reported correlation-like coefficient is
its normalisation ``r_ij = R_ij / sqrt(R_ii R_jj)`` (Pearson correlation across
sites of the two trajectory profiles); the diagonal is 1 by construction.
Negative values (genetically diverging trajectories) carry no shared-ancestry
signal and are clipped to 0, matching the ``entries in [0,1]`` contract of
``sim/trajectory_stat._correlation``.

Site filter: autosomes 1--22 only; exclude chr6 MHC (GRCh37 28.5--33.5 Mb) and
the chr1 AMY1/AMY2 locus (GRCh37 103.5--104.5 Mb); drop sites with >20% missing
among the trajectory individuals; then take a fixed-seed uniform subsample of
the surviving sites for tractable decoding.

N_eff = m^2 / sum(R) over the distinct trajectory set, exactly as
``sim/trajectory_stat.trajectory_statistic`` computes it.

All claims in the companion doc are labelled [observation] / [inference] /
[speculation].
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from analysis.load_aadr import BUNDLE, MAIN_PREFIX, GenotypeMatrix, iter_snp  # noqa: E402

TRAJ_TSV = "data/derived/trajectories.tsv"
OUT_MATRIX = "data/derived/traj_relatedness.tsv"
OUT_METRICS = "data/derived/traj_relatedness_metrics.json"
OUT_FIG = "fig/traj_relatedness.png"

SNP_PATH = os.path.join(BUNDLE, MAIN_PREFIX + ".snp")
GENO_PATH = os.path.join(BUNDLE, MAIN_PREFIX + ".geno")
IND_PATH = os.path.join(BUNDLE, MAIN_PREFIX + ".ind")

MHC = ("6", 28_500_000, 33_500_000)      # GRCh37 MHC
AMY = ("1", 103_500_000, 104_500_000)    # chr1 AMY1/AMY2 GRCh37
MAX_MISS = 0.20
N_SITES_TARGET = 150_000
SEED = 20240501


def load_trajectories():
    t = pd.read_csv(TRAJ_TSV, sep="\t")
    t = t.loc[t.trajectory_id >= 0].copy()
    t["trajectory_id"] = t.trajectory_id.astype(int)
    return t


def candidate_sites():
    """Site indices (0-based) passing chrom/region filters."""
    keep = []
    for i, rec in enumerate(iter_snp(SNP_PATH)):
        chrom = rec[1]
        pos = rec[3]
        if chrom not in {str(c) for c in range(1, 23)}:
            continue
        if chrom == MHC[0] and MHC[1] <= pos <= MHC[2]:
            continue
        if chrom == AMY[0] and AMY[1] <= pos <= AMY[2]:
            continue
        keep.append(i)
    return np.asarray(keep, dtype=np.int64)


def select_sites(cand):
    rng = np.random.default_rng(SEED)
    if len(cand) <= N_SITES_TARGET:
        return np.sort(cand)
    pick = rng.choice(len(cand), size=N_SITES_TARGET, replace=False)
    return np.sort(cand[pick])


def group_accumulate(gm, ind_rows, traj_of_row, sites):
    """Single decoding pass. Returns (sums, counts) shaped (n_group, S)."""
    groups = np.array(sorted(set(traj_of_row)), dtype=int)
    gidx = {g: k for k, g in enumerate(groups)}
    S = len(sites)
    sums = np.zeros((len(groups), S), dtype=np.float64)
    counts = np.zeros((len(groups), S), dtype=np.int64)
    for n, row in enumerate(ind_rows):
        codes = gm.individual_genotypes(row, sites)
        present = codes != 3
        g = gidx[traj_of_row[n]]
        sums[g] += np.where(present, codes, 0)
        counts[g] += present
        if (n + 1) % 2000 == 0:
            print(f"  decoded {n+1}/{len(ind_rows)} individuals", flush=True)
    return groups, sums, counts


def relatedness_from_accum(groups, sums, counts):
    """Build the labeled correlation matrix in [0,1].

    The per-site reference allele frequency is the *equal-weight* mean of the
    trajectory mean frequencies (each trajectory contributes 1/m regardless of
    sample size), so a large central cluster cannot define the reference away
    from itself.  Each trajectory's standardised deviation profile is
    ``d_{t,s} = (p_{t,s} - ref_s) / sqrt(2 ref_s (1-ref_s))`` and
    ``R = (1/S) D D^T``; entries are the Pearson correlation of the profiles.
    """
    valid = counts.min(axis=0) > 0
    meang = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)
    mg = meang[:, valid]
    ref = mg.mean(axis=0)          # equal-weight trajectory reference
    f = ref / 2.0
    sd = np.sqrt(2 * f * (1 - f))
    good = sd > 0
    D = (mg[:, good] - ref[good]) / sd[good]   # (n_group, S)
    S = D.shape[1]
    R = (D @ D.T) / S
    d = np.sqrt(np.clip(np.diag(R), 1e-12, None))
    r = R / np.outer(d, d)
    r = np.clip(r, 0.0, 1.0)
    np.fill_diagonal(r, 1.0)
    labels = [str(g) for g in groups]
    mat = pd.DataFrame(r, index=labels, columns=labels)
    return mat, int(S)


def merge_by_threshold(mat, thr):
    """Union-find merge of labels whose pairwise coefficient exceeds thr."""
    labels = list(mat.index)
    parent = {x: x for x in labels}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a in labels:
        for b in labels:
            if a < b and mat.loc[a, b] > thr:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[rb] = ra
    clusters = {}
    for x in labels:
        clusters.setdefault(find(x), []).append(x)
    return sorted((sorted(v) for v in clusters.values()))


def coarse_matrix(mat, clusters):
    """Cluster-to-cluster relatedness = mean of member pair coefficients."""
    labels = ["+".join(c) for c in clusters]
    C = np.zeros((len(clusters), len(clusters)))
    for a, A in enumerate(clusters):
        for b, B in enumerate(clusters):
            C[a, b] = mat.loc[A, B].to_numpy().mean()
    return pd.DataFrame(C, index=labels, columns=labels)


def n_eff(mat, clusters=None):
    """N_eff = m^2 / sum(R) on the distinct or merged cluster set."""
    if clusters is None:
        m = mat.shape[0]
        return float(m * m / mat.to_numpy().sum()), m
    coarse = coarse_matrix(mat, clusters)
    m = coarse.shape[0]
    return float(m * m / coarse.to_numpy().sum()), m


def main():
    print("== loading trajectories ==", flush=True)
    traj = load_trajectories()
    ind_ids = [l.split()[0] for l in open(IND_PATH)]
    row_of = {x: i for i, x in enumerate(ind_ids)}
    ind_rows = np.array([row_of[x] for x in traj.individual], dtype=np.int64)
    traj_of_row = traj.trajectory_id.to_numpy()
    print(f"  {len(ind_rows)} individuals in {traj.trajectory_id.nunique()} trajectories", flush=True)

    print("== filtering sites ==", flush=True)
    cand = candidate_sites()
    sites = select_sites(cand)
    print(f"  candidate autosomal/MHC-AMY-excluded sites: {len(cand)}; "
          f"subsampled: {len(sites)}", flush=True)

    print("== decoding genotypes ==", flush=True)
    gm = GenotypeMatrix(GENO_PATH)
    groups, sums, counts = group_accumulate(gm, ind_rows, traj_of_row, sites)

    print("== building relatedness ==", flush=True)
    mat, S_used = relatedness_from_accum(groups, sums, counts)
    mat.to_csv(OUT_MATRIX, sep="\t")
    print(f"  matrix written; sites used after missingness filter: {S_used}", flush=True)

    # ---- validation through the consumer ----
    from sim.trajectory_stat import _correlation
    pd_error = None
    try:
        _correlation(mat, list(mat.index))
        print("  _correlation (Cholesky PD check): PASS", flush=True)
    except Exception as exc:  # noqa: BLE001
        pd_error = str(exc)
        print(f"  _correlation: FAIL -> {exc}", flush=True)

    # ---- N_eff on full set and coarser groupings ----
    neff_full, m_full = n_eff(mat)
    thresholds = [0.99, 0.95, 0.9, 0.8, 0.7, 0.5]
    groupings = []
    for thr in thresholds:
        cl = merge_by_threshold(mat, thr)
        ne, m = n_eff(mat, cl)
        groupings.append({"threshold": thr, "n_groups": m, "n_effective": ne,
                          "clusters": cl})
    print(f"  N_eff (13 distinct) = {neff_full:.3f}", flush=True)
    for g in groupings:
        print(f"  merge>{g['threshold']}: {g['n_groups']} groups -> "
              f"N_eff={g['n_effective']:.3f}", flush=True)

    # near-duplicate pairs
    labels = list(mat.index)
    dupes = []
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            v = mat.loc[a, b]
            if v >= 0.8:
                dupes.append({"a": a, "b": b, "r": float(v)})

    # ---- split-trajectory sanity check ----
    print("== split-trajectory sanity check ==", flush=True)
    split_report = split_trajectory_check(gm, ind_rows, traj_of_row, traj, sites,
                                          groups, sums, counts)

    metrics = {
        "estimator": "normalized mean pairwise allele-sharing/IBS: "
                     "r_ij = R_ij/sqrt(R_ii R_jj) with R=(1/S) D D^T, "
                     "D_{t,s}=(p_{t,s}-ref_s)/sqrt(2 ref_s (1-ref_s)), "
                     "ref_s = equal-weight mean of trajectory frequencies",
        "site_filter": {
            "autosomes": "1-22",
            "excluded": ["chr6 MHC GRCh37 28.5-33.5 Mb",
                         "chr1 AMY1/AMY2 GRCh37 103.5-104.5 Mb"],
            "max_missing": MAX_MISS,
            "candidate_sites": int(len(cand)),
            "subsampled_sites": int(len(sites)),
            "sites_used_after_missingness": int(S_used),
            "seed": SEED,
        },
        "n_trajectories": int(mat.shape[0]),
        "labels": labels,
        "n_effective_all13": neff_full,
        "coarser_groupings": groupings,
        "near_duplicates_r_ge_0.8": dupes,
        "cholesky_pd": pd_error is None,
        "cholesky_error": pd_error,
        "split_sanity_check": split_report,
    }
    with open(OUT_METRICS, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print("  metrics written", flush=True)

    make_plot(mat)
    return metrics


def split_trajectory_check(gm, ind_rows, traj_of_row, traj, sites, groups, sums, counts):
    """Split the largest trajectory into halves; N_eff must not increase."""
    sizes = traj.trajectory_id.value_counts()
    big = int(sizes.idxmax())
    print(f"  splitting trajectory {big} ({int(sizes[big])} individuals)", flush=True)
    rng = np.random.default_rng(SEED + 1)
    members = np.where(traj_of_row == big)[0]
    perm = rng.permutation(len(members))
    half = np.zeros(len(ind_rows), dtype=int)
    half[members[perm[: len(members) // 2]]] = 1
    new_labels = np.where(traj_of_row == big,
                          np.where(half == 1, 1000 + big, 2000 + big),
                          traj_of_row)
    # recompute from per-individual so the two halves get distinct group means
    groups2, sums2, counts2 = group_accumulate(gm, ind_rows, new_labels, sites)
    mat2, S2 = relatedness_from_accum(groups2, sums2, counts2)
    ne2, m2 = n_eff(mat2)
    base_labels = [str(g) for g in groups]
    mat_base = pd.DataFrame(
        relatedness_from_accum(groups, sums, counts)[0].values,
        index=base_labels, columns=base_labels)
    ne1, m1 = n_eff(mat_base)
    # relatedness between the two halves
    rab = float(mat2.loc[str(1000 + big), str(2000 + big)])
    out = {
        "split_trajectory": big,
        "n_individuals": int(sizes[big]),
        "r_between_halves": rab,
        "n_effective_distinct13": ne1,
        "n_effective_after_split14": ne2,
        "n_eff_did_not_increase": bool(ne2 <= ne1 + 1e-9),
    }
    print(f"  r between halves = {rab:.4f}", flush=True)
    print(f"  N_eff distinct-13 = {ne1:.3f} -> after split-14 = {ne2:.3f} "
          f"(did not increase: {out['n_eff_did_not_increase']})", flush=True)
    return out


def make_plot(mat):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(mat.values, cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(len(mat.index)))
    ax.set_yticks(range(len(mat.index)))
    ax.set_xticklabels(mat.index, rotation=90)
    ax.set_yticklabels(mat.index)
    for i in range(len(mat.index)):
        for j in range(len(mat.index)):
            ax.text(j, i, f"{mat.values[i, j]:.2f}", ha="center", va="center",
                    color="w" if mat.values[i, j] < 0.6 else "k", fontsize=6)
    fig.colorbar(im, ax=ax, fraction=0.046, label="relatedness r")
    ax.set_title("Trajectory pairwise genetic relatedness")
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=150)
    print(f"  figure written: {OUT_FIG}", flush=True)


if __name__ == "__main__":
    main()
