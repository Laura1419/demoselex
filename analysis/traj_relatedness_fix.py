#!/usr/bin/env python3
"""Trajectory relatedness R and effective independent trajectory count N_eff.

Estimator (documented): for each autosomal site s (excluding MHC and AMY1/AMY2),
compute each trajectory's mean allele frequency f_{k,s} from the called AADR 1240K
genotypes of its members.  Remove the shared site allele-frequency spectrum
(center across trajectories) and take the correlation across sites between
trajectory frequency profiles:

    C_ij = corr_s( f_i. , f_j. )

This is the genome-wide correlation of trajectory allele-frequency *deviations*,
i.e. the shared-drift covariance that the GLS meta-analysis needs.  Correlation
entries may be negative (anti-correlated profiles); the consumer
(sim/trajectory_stat.py::_correlation) requires a symmetric matrix with unit
diagonal and entries in [0,1] that is positive definite, so we clip negatives to
zero and project to the nearest positive-definite correlation matrix by
eigenvalue flooring, then rescale to unit diagonal.  Every transformation is
reported so the choice is auditable.

N_eff = m^2 / sum(R) over the distinct trajectory set.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from analysis.load_aadr import GenotypeMatrix, iter_snp, parse_ind  # noqa: E402

ROOT = Path("/home/erikg/demoselex")
BUNDLE = ROOT / "data/raw/aadr-v66.p1"
GENO = str(BUNDLE / "v66.p1_1240K.aadr.patch.PUB.geno")
SNP = str(BUNDLE / "v66.p1_1240K.aadr.patch.PUB.snp")
IND = str(BUNDLE / "v66.p1_1240K.aadr.patch.PUB.ind")
TRAJ = ROOT / "data/derived/trajectories.tsv"
OUT_MATRIX = ROOT / "data/derived/traj_relatedness.tsv"
OUT_METRICS = ROOT / "data/derived/traj_relatedness_metrics.json"

N_SITES = 150_000
SEED = 20240501
MAX_MISSING = 0.5


def site_index() -> np.ndarray:
    chrom, pos = [], []
    for rec in iter_snp(SNP):
        chrom.append(rec.chrom if hasattr(rec, "chrom") else rec[1])
        pos.append(rec.pos if hasattr(rec, "pos") else rec[3])
    chrom = np.asarray(chrom, dtype=object)
    pos = np.asarray(pos, dtype=np.int64)
    ok = np.isin(chrom, [str(c) for c in range(1, 23)])
    mhc = (chrom == "6") & (pos >= 28_500_000) & (pos <= 33_500_000)
    amy = (chrom == "1") & (pos >= 103_500_000) & (pos <= 104_500_000)
    ok &= ~(mhc | amy)
    cand = np.flatnonzero(ok)
    rng = np.random.default_rng(SEED)
    take = min(N_SITES, cand.size)
    return np.sort(rng.choice(cand, size=take, replace=False))


def main() -> int:
    traj = pd.read_csv(TRAJ, sep="\t")
    # trajectory_id == -1 is the "unassigned" bucket, not a trajectory
    traj = traj[traj.trajectory_id >= 0]
    ind_ids = [r.id for r in parse_ind(IND)]
    row_of = {x: i for i, x in enumerate(ind_ids)}
    traj = traj[traj.individual.isin(row_of)].reset_index(drop=True)
    groups = sorted(traj.trajectory_id.unique().tolist())
    gidx = {g: k for k, g in enumerate(groups)}
    rows_by_group: dict[int, list[int]] = {g: [] for g in groups}
    for ind, g in zip(traj.individual, traj.trajectory_id):
        rows_by_group[g].append(row_of[ind])

    sites = site_index()
    gm = GenotypeMatrix(GENO)
    n_g = len(groups)
    n_s = sites.size
    sums = np.zeros((n_g, n_s), dtype=np.float32)
    counts = np.zeros((n_g, n_s), dtype=np.float32)

    for g in groups:
        k = gidx[g]
        for i in rows_by_group[g]:
            codes = gm.individual_genotypes(i, sites)
            called = codes <= 2
            sums[k] += np.where(called, codes, 0)
            counts[k] += called
        print(f"  trajectory {g}: {len(rows_by_group[g])} individuals", flush=True)

    with np.errstate(invalid="ignore", divide="ignore"):
        freq = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    # keep sites where every trajectory has a frequency so the profile is defined
    complete = np.isfinite(freq).all(axis=0)
    F = freq[:, complete]
    S = F.shape[1]
    print(f"  sites usable (all trajectories covered): {S}", flush=True)

    # center each site across trajectories -> remove the shared MAF spectrum
    dev = F - F.mean(axis=0, keepdims=True)
    # correlation of profiles across sites
    with np.errstate(invalid="ignore", divide="ignore"):
        C = (dev @ dev.T) / S
        d = np.sqrt(np.clip(np.diag(C), 1e-12, None))
        C = C / np.outer(d, d)
    C = np.clip(C, -1.0, 1.0)

    # --- clip negatives, project to nearest PD correlation, rescale ---
    C0 = np.clip(C, 0.0, 1.0)
    np.fill_diagonal(C0, 1.0)
    w, V = np.linalg.eigh(C0)
    w_clipped = np.clip(w, 1e-6, None)
    Cpd = V @ np.diag(w_clipped) @ V.T
    dd = np.sqrt(np.clip(np.diag(Cpd), 1e-12, None))
    Cpd = Cpd / np.outer(dd, dd)
    Cpd = np.clip(Cpd, 0.0, 1.0)
    np.fill_diagonal(Cpd, 1.0)
    # symmetric
    Cpd = (Cpd + Cpd.T) / 2.0

    labels = [str(g) for g in groups]
    mat = pd.DataFrame(Cpd, index=labels, columns=labels)
    mat.to_csv(OUT_MATRIX, sep="\t")

    # validate with the consumer
    from sim.trajectory_stat import _correlation
    pd_ok, pd_err = True, None
    try:
        _correlation(mat, labels)
    except Exception as exc:  # noqa: BLE001
        pd_ok, pd_err = False, str(exc)

    m = len(labels)
    n_eff = float(m * m / mat.to_numpy().sum())

    # split sanity check: split the largest trajectory in half -> N_eff must not rise
    big = int(traj.trajectory_id.value_counts().idxmax())
    members = rows_by_group[big]
    half = len(members) // 2
    split_sums = np.zeros((n_g + 1, n_s), dtype=np.float32)
    split_counts = np.zeros((n_g + 1, n_s), dtype=np.float32)
    for k, g in enumerate(groups):
        for i in rows_by_group[g]:
            codes = gm.individual_genotypes(i, sites)
            called = codes <= 2
            split_sums[k] += np.where(called, codes, 0)
            split_counts[k] += called

    def profile(counts_m, sums_m):
        with np.errstate(invalid="ignore", divide="ignore"):
            f = np.divide(sums_m, counts_m, out=np.full_like(sums_m, np.nan),
                          where=counts_m > 0)
        cm = np.isfinite(f).all(axis=0)
        Fm = f[:, cm]
        dm = Fm - Fm.mean(axis=0, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            Cm = (dm @ dm.T) / Fm.shape[1]
            dg = np.sqrt(np.clip(np.diag(Cm), 1e-12, None))
            Cm = Cm / np.outer(dg, dg)
        Cm = np.clip((Cm + Cm.T) / 2, 0.0, 1.0)
        np.fill_diagonal(Cm, 1.0)
        return Cm

    half_sums = split_sums.copy()
    half_counts = split_counts.copy()
    # replace the big group's row with half of its individuals, add a new row for the rest
    a = np.zeros(n_s, dtype=np.float32); b = np.zeros(n_s, dtype=np.float32)
    ca = np.zeros(n_s, dtype=np.float32); cb = np.zeros(n_s, dtype=np.float32)
    for j, i in enumerate(members):
        codes = gm.individual_genotypes(i, sites)
        called = codes <= 2
        tgt_s, tgt_c = (a, ca) if j < half else (b, cb)
        tgt_s += np.where(called, codes, 0); tgt_c += called
    new_sums = np.vstack([half_sums, b[None, :]])
    new_counts = np.vstack([half_counts, cb[None, :]])
    new_sums[gidx[big]] = a
    new_counts[gidx[big]] = ca
    Csplit = profile(new_counts, new_sums)
    m2 = Csplit.shape[0]
    n_eff_split = float(m2 * m2 / Csplit.sum())

    metrics = {
        "estimator": "correlation of trajectory allele-frequency deviations across "
                     f"{S} autosomal sites (MHC and AMY excluded); negatives clipped to 0, "
                     "eigenvalue-floored to the nearest PD correlation, unit diagonal rescaled",
        "n_trajectories": m,
        "n_individuals": int(len(traj)),
        "sites_used": int(S),
        "seed": SEED,
        "raw_correlation_min": float(C.min()),
        "raw_correlation_max_offdiag": float((C - np.eye(m)).max()),
        "n_effective_all": n_eff,
        "cholesky_pd": pd_ok,
        "cholesky_error": pd_err,
        "split_sanity_check": {
            "split_trajectory": big,
            "n_effective_before": n_eff,
            "n_effective_after_split": n_eff_split,
            "n_eff_did_not_increase": bool(n_eff_split <= n_eff + 1e-9),
        },
        "pairwise_ge_0.5": [
            {"a": labels[i], "b": labels[j], "r": float(mat.iloc[i, j])}
            for i in range(m) for j in range(i + 1, m) if mat.iloc[i, j] >= 0.5
        ],
    }
    OUT_METRICS.write_text(json.dumps(metrics, indent=1))

    print(f"  N_eff(all {m}) = {n_eff:.3f}  PD={pd_ok}")
    print("  matrix:")
    print(mat.round(3).to_string())
    print(f"  split check: N_eff {n_eff:.3f} -> {n_eff_split:.3f} "
          f"(did not increase: {metrics['split_sanity_check']['n_eff_did_not_increase']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
