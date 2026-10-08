#!/usr/bin/env python3
"""Between-history ordered-rise test (pre-registered in docs/traj-trend-prereg.md).

Tests a THIRD regime distinct from within-trajectory change: each history is
internally flat, but successive histories (in time order) carry higher mean allele
frequency. Statistic = relatedness-aware GLS slope of history mean frequency on
history date, at the level of the merged replication groups.

Pre-registration: docs/traj-trend-prereg.md (written before this ran).
Outputs: data/derived/traj_trend.tsv (written incrementally per genotype chunk).

Usage:  .venv/bin/python analysis/traj_trend.py --max-variants 2000
"""

from __future__ import annotations

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "analysis"))
sys.path.insert(0, str(REPO / "sim"))

import load_akbari_geno as akg  # noqa: E402
from tier2_scan import build_sample_map  # noqa: E402
from trajectory_stat import _correlation  # noqa: E402

SIG_ABS_X = 3.16
SEED = 42
POOL_CACHE = REPO / "data/derived/tier2_pool_cache.pkl"
REL_V3 = REPO / "data/derived/traj_relatedness_v3.tsv"
TRAJ_V3 = REPO / "data/derived/trajectories_v3.tsv"

COLUMNS = [
    "locus", "chrom", "pos", "ref", "alt", "akbari_X", "akbari_AF",
    "missingness", "m", "beta", "se", "z", "ordering_concordant",
    "ordering_possible", "relatedness_frac", "p_perm", "n_perm",
    "is_candidate", "source", "n_controls", "p_matched", "bh_fdr",
]


# --------------------------------------------------------------------------- #
# Panel: reuse the prepared QC-confirmed cache (copied into data/derived).
# --------------------------------------------------------------------------- #
def load_pool():
    return pickle.loads(POOL_CACHE.read_bytes())


def select_panel(pool, max_variants):
    pool = pool[pool["pgen_index"] >= 0].reset_index(drop=True)
    rng = np.random.default_rng(SEED)
    cand_idx = np.nonzero(pool["is_candidate"].to_numpy())[0]
    bg_idx = np.nonzero(~pool["is_candidate"].to_numpy())[0]
    n_cand = min(len(cand_idx), max(1, max_variants // 4))
    if len(cand_idx) > n_cand:
        cand_idx = np.sort(rng.choice(cand_idx, size=n_cand, replace=False))
    n_bg = min(len(bg_idx), max(0, max_variants - len(cand_idx)))
    if len(bg_idx) > n_bg:
        bg_idx = np.sort(rng.choice(bg_idx, size=n_bg, replace=False))
    panel = pool.iloc[np.concatenate([cand_idx, bg_idx])].reset_index(drop=True)
    return panel.sort_values("pgen_index").reset_index(drop=True)


# --------------------------------------------------------------------------- #
# GLS machinery (fixed design -> linear functional of p)
# --------------------------------------------------------------------------- #
def gls_design(r_mat):
    """Return c, denom, tc for beta_hat = c'p / denom (t from -date_BP)."""
    rinv = np.linalg.inv(r_mat)
    return rinv


def gls_z_chunk(p, tc, rinv, denom):
    """Vectorised GLS z for p of shape (n, m); tc (m,), rinv (m,m), denom scalar."""
    c = rinv @ tc
    beta = (p @ c) / denom
    ones_r = rinv.sum(axis=1)  # R^{-1} 1
    mu = (p @ ones_r) / (ones_r.sum())
    pc = p - mu[:, None]
    r = pc - beta[:, None] * tc[None, :]
    rinv_r = r @ rinv.T
    sig2 = np.einsum("ij,ij->i", r, rinv_r) / (p.shape[1] - 2)
    se = np.sqrt(sig2 / denom)
    with np.errstate(divide="ignore", invalid="ignore"):
        z = beta / se
    return beta, se, z


def ordering_concordant(p):
    """# successive increases in p when histories sorted by ascending t."""
    order = np.argsort(T_GLOBAL)  # fixed group order by date
    ps = p[:, order]
    return (np.diff(ps, axis=1) > 0).sum(axis=1)


def permutation_p(p, beta_obs, se_obs, tc, rinv, denom):
    """Exchangeable permutation P holding the observed scale ``se_obs`` fixed.

    Re-estimating sigma^2 per permutation is not exchangeable (sigma^2 is
    minimised at the observed arrangement), so the permutation distribution must
    use a common scale; this is the fixed-scale permutation null of the prereg.
    """
    m = p.shape[1]
    import itertools
    perms = list(itertools.permutations(range(m)))
    c = rinv @ tc
    az = np.abs(beta_obs / se_obs)
    counts = np.zeros(p.shape[0], dtype=np.int64)
    for perm in perms:
        bp = (p[:, perm] @ c) / denom
        counts += (np.abs(bp / se_obs) >= az).astype(np.int64)
    return (1 + counts) / (1 + len(perms)), len(perms)


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return out
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(m)
    tmp[order] = q
    out[ok] = tmp
    return out


# --------------------------------------------------------------------------- #
# Matched-control empirical null (same matcher family as tier2_scan)
# --------------------------------------------------------------------------- #
def matched_null(df, query_is_candidate, min_controls=20):
    """Matched-control empirical P for every query locus with finite z.

    Controls are always the background (``is_candidate == False``) loci; queries
    are the loci selected by ``query_is_candidate`` (None = all finite loci).
    """
    bg = df[(~df["is_candidate"]) & np.isfinite(df["z"])]
    qmask = (np.isfinite(df["z"])
             & (query_is_candidate if query_is_candidate is not None
                else pd.Series(True, index=df.index)))
    cand = df[qmask]
    p_matched = np.full(len(cand), np.nan)
    n_ctrl = np.zeros(len(cand), dtype=int)
    for j, (_i, crow) in enumerate(cand.iterrows()):
        pool = bg[(bg["chrom"] == crow["chrom"])
                  & ((bg["pos"] - crow["pos"]).abs() >= 1_000_000)]
        chosen = None
        for scale in (1.0, 2.0, 4.0):
            sel = pool[((pool["akbari_AF"] - crow["akbari_AF"]).abs() <= min(0.10 * scale, 0.5))
                       & ((pool["missingness"] - crow["missingness"]).abs() <= min(0.10 * scale, 0.5))]
            if len(sel) >= min_controls:
                chosen = sel
                break
        if chosen is None:
            chosen = pool
        n_ctrl[j] = len(chosen)
        if len(chosen) == 0:
            continue
        zs = np.abs(chosen["z"].to_numpy(float))
        p_matched[j] = (1 + int((zs >= abs(crow["z"])).sum())) / (len(chosen) + 1)
    return cand.index, p_matched, n_ctrl


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--max-variants", type=int, default=8000)
    ap.add_argument("--chunk-variants", type=int, default=1000)
    ap.add_argument("--min-controls", type=int, default=20)
    ap.add_argument("--out-tsv", type=Path,
                    default=REPO / "data/derived/traj_trend.tsv")
    args = ap.parse_args(argv)

    t0 = time.time()
    smap = build_sample_map(verbose=False)
    rel = pd.read_csv(REL_V3, sep="\t", index_col=0)
    group_order = list(rel.index)
    r_mat = _correlation(rel, group_order)
    m = len(group_order)

    # Fixed temporal placement: t = -mean(date_BP) per group (larger = more recent).
    gsizes = smap.groupby("merged_group")["date_BP"].mean()
    t_chron = -np.array([gsizes[g] for g in group_order], float)
    t_c = t_chron - t_chron.mean()
    rinv = gls_design(r_mat)
    denom = float(t_c @ rinv @ t_c)
    rel_frac = 1.0 - denom / float(t_c @ t_c)
    global T_GLOBAL
    T_GLOBAL = t_chron

    print(f"units: m={m} merged groups {group_order}; N_eff={m*m/r_mat.sum():.3f}; "
          f"t_chron={np.round(t_chron,0).tolist()}; relatedness_frac={rel_frac:.3f}",
          flush=True)

    pool = load_pool()
    panel = select_panel(pool, args.max_variants)
    n_panel = len(panel)
    print(f"panel: {n_panel} variants "
          f"(candidates={int(panel['is_candidate'].sum())}, "
          f"background={int((~panel['is_candidate']).sum())})", flush=True)

    subset = smap["pgen_full_index"].to_numpy(np.int32)
    n_ind = len(subset)
    gord = {g: i for i, g in enumerate(group_order)}
    gid = smap["merged_group"].map(gord).to_numpy(np.int64)
    sorter = np.argsort(gid, kind="stable")
    gr_sorted = gid[sorter]
    starts = np.searchsorted(gr_sorted, np.arange(m))

    labels = (panel["CHROM"].astype(str) + ":" + panel["POS"].astype(str) + ":"
              + panel["REF"].astype(str) + ":" + panel["ALT"].astype(str)).to_numpy()

    out = args.out_tsv
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\t".join(COLUMNS) + "\n")

    geno = akg.AkbariGeno(sample_subset=subset)
    assert geno.n_samples == n_ind

    for start in range(0, n_panel, args.chunk_variants):
        stop = min(start + args.chunk_variants, n_panel)
        pgen_idx = panel["pgen_index"].to_numpy()[start:stop]
        dos = geno.read_variants(pgen_idx)
        finite = np.isfinite(dos)
        d0 = np.where(finite, dos, 0.0)
        gsum = np.add.reduceat(d0[:, sorter], starts, axis=1).T       # (m, chunk)
        gcnt = np.add.reduceat(finite[:, sorter], starts, axis=1).T   # (m, chunk)
        with np.errstate(invalid="ignore"):
            p = np.where(gcnt > 0, gsum / (2.0 * gcnt), np.nan)  # (m, chunk) AF
        miss = 1.0 - (finite.sum(axis=1) / float(n_ind))
        # allele frequency matrix as (chunk, m)
        pm = p.T
        ok = np.isfinite(pm).all(axis=1)

        beta = np.full(stop - start, np.nan)
        se = np.full(stop - start, np.nan)
        z = np.full(stop - start, np.nan)
        ordc = np.full(stop - start, np.nan)
        pp = np.full(stop - start, np.nan)
        nperm = 0
        if ok.any():
            sub = pm[ok]
            b, s, zz = gls_z_chunk(sub, t_c, rinv, denom)
            beta[ok] = b
            se[ok] = s
            z[ok] = zz
            ordc[ok] = ordering_concordant(sub)
            ppv, nperm = permutation_p(sub, b, s, t_c, rinv, denom)
            pp[ok] = ppv

        chunk_out = pd.DataFrame({
            "locus": labels[start:stop],
            "chrom": panel["CHROM"].to_numpy()[start:stop],
            "pos": panel["POS"].to_numpy()[start:stop],
            "ref": panel["REF"].to_numpy()[start:stop],
            "alt": panel["ALT"].to_numpy()[start:stop],
            "akbari_X": panel["X"].to_numpy()[start:stop],
            "akbari_AF": panel["AF"].to_numpy()[start:stop],
            "missingness": miss,
            "m": np.full(stop - start, m),
            "beta": beta, "se": se, "z": z,
            "ordering_concordant": ordc,
            "ordering_possible": np.full(stop - start, m - 1),
            "relatedness_frac": np.full(stop - start, rel_frac),
            "p_perm": pp, "n_perm": np.full(stop - start, nperm),
            "is_candidate": panel["is_candidate"].to_numpy()[start:stop],
            "source": np.where(panel["is_candidate"].to_numpy()[start:stop],
                               "akbari_sig", "background"),
            "n_controls": np.nan, "p_matched": np.nan, "bh_fdr": np.nan,
        })
        chunk_out.to_csv(out, sep="\t", index=False, header=False, mode="a")
        print(f"  chunk {start}-{stop}: {int(ok.sum())} scored, "
              f"t={time.time()-t0:.0f}s", flush=True)

    df = pd.read_csv(out, sep="\t")
    cond = (df["missingness"] <= 0.25) & np.isfinite(df["p_perm"])
    # Calibrated null = matched-control empirical P (prereg decision rule + gate).
    idx, pm_, nc = matched_null(df, cond, args.min_controls)
    df.loc[idx, "p_matched"] = pm_
    df.loc[idx, "n_controls"] = nc
    df.loc[cond, "bh_fdr"] = bh_fdr(df.loc[cond, "p_matched"])
    df.to_csv(out, sep="\t", index=False)

    # ---- summary --------------------------------------------------------- #
    d = df[cond].copy()
    d["absX"] = d["akbari_X"].abs()
    cal = float((d.loc[~d["is_candidate"], "p_matched"] <= 0.05).mean())
    lo = d[~d["is_candidate"]]
    hi = d[d["is_candidate"]]
    print(f"scored={len(d)}; calibration: background p_matched<=0.05 frac={cal:.3f}")
    print(f"low|X| frac p_matched<=0.05 = {float((lo['p_matched']<=0.05).mean()):.3f} "
          f"(n={len(lo)}); high|X| frac = {float((hi['p_matched']<=0.05).mean()):.3f} "
          f"(n={len(hi)})")
    print(f"low|X| frac p_perm<=0.05 = {float((lo['p_perm']<=0.05).mean()):.3f}; "
          f"high|X| = {float((hi['p_perm']<=0.05).mean()):.3f}")
    top = d.nsmallest(10, "p_matched")[
        ["locus", "akbari_X", "beta", "z", "ordering_concordant", "p_matched", "p_perm"]]
    print(top.to_string(index=False))
    print(f"done t={time.time()-t0:.0f}s -> {out}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
