#!/usr/bin/env python3
"""Stage A v2: ordered-path temporal trajectories from genome-wide similarity + time.

Supersedes ``analysis/trajectory_infer.py`` (v1).  v1 (time-window k-means) failed:
one blob of 9,648/10,369 individuals spanning 0-65,000 BP plus fragments of size
5 and 2 (see ``docs/traj-validation.md``).  v2 makes three structural changes
required by that failure:

1. Modern individuals (``date_method`` starts with "Modern", or ``date_bp == 0``)
   never enter trajectory inference.  Only non-modern individuals with a usable
   date, usable coordinates, and a sane date window are used.
2. A per-individual coverage floor (fraction of selected sites called) removes
   samples too sparse to place reliably.
3. Trajectories are **ordered paths**, not clusters.  The time axis is cut into
   overlapping bins; each bin is clustered (seeded k-means) in PCA space; clusters
   are linked across ADJACENT bins by genetic centroid distance under a hard gate
   with a greedy one-to-one assignment; a trajectory is a chain of linked
   clusters.  Chains are capped (min size, max span) and no chain skips a bin.

The embedding builder accepts an explicit exclusion list of chromosome regions;
the default excludes chr6 MHC and chr1 AMY1/AMY2 (GRCh37).

Genotype I/O uses the AADR v66 TGENO layout (48-byte header, individual-major,
2 bits/site, 4 sites/byte, code 3 = missing).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time

import numpy as np

BUNDLE = "data/raw/aadr-v66.p1"
PREFIX = "v66.p1_1240K.aadr.patch.PUB"
ANNO = "data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.anno"
GENO = f"{BUNDLE}/{PREFIX}.geno"
IND = f"{BUNDLE}/{PREFIX}.ind"
SNP = f"{BUNDLE}/{PREFIX}.snp"

# GRCh37. Default exclusion: chr6 MHC and chr1 AMY1/AMY2.
DEFAULT_EXCLUDE = (("6", 28_500_000, 33_500_000), ("1", 103_500_000, 104_500_000))

# Known-structure annotation keywords (post-hoc only; never inputs).
WHG_KEYS = ("whg", "mesolithic", "loschbour", "villabruna", "kozlowski",
            "hunter", "motala", "satsurblia", "karelia", "iron_gates",
            "serbia_meso", "bichon", "oberkassel", "la_bra", "el_mena")
EEF_KEYS = ("eef", "early_farmer", "earlyeuropeanfarmer", "lbk", "anatolia",
            "barcin", "tepecik", "starcevo", "cardial", "impressa", "greece_en",
            "boncuklu", "aegean", "karsdorf", "neolithic")
STEPPE_KEYS = ("steppe", "yamnaya", "corded", "afanasievo", "andronovo",
               "srubnaya", "poltavka", "catacomb", "sintashta", "bell_beaker",
               "baden", "vucedol", "cw_", "globular")


# --------------------------------------------------------------------------- #
# metadata
# --------------------------------------------------------------------------- #
def read_ind(path=IND):
    ids = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if p:
                ids.append(p[0])
    return ids


def read_anno(path=ANNO):
    """Genetic ID -> dict(date, method, lat, lon, group_id, locality)."""
    info = {}
    with open(path, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr)
        for row in rdr:
            gid = row[0]
            if gid in info:
                continue
            def fl(i):
                v = row[i].strip()
                try:
                    return float(v)
                except (ValueError, IndexError):
                    return None
            info[gid] = {
                "date": fl(10), "method": row[9], "lat": fl(17), "lon": fl(18),
                "group_id": row[14], "locality": row[15],
            }
    return info


def read_snp(path=SNP):
    chrom, pos = [], []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            chrom.append(p[1])
            pos.append(int(p[3]))
    return np.array(chrom), np.array(pos, dtype=np.int64)


def excluded_mask(chrom, pos, regions):
    bad = np.zeros(len(pos), dtype=bool)
    for c, lo, hi in regions:
        bad |= (chrom == str(c)) & (pos >= lo) & (pos < hi)
    return bad


# --------------------------------------------------------------------------- #
# genotype I/O (TGENO)
# --------------------------------------------------------------------------- #
class TGENO:
    def __init__(self, path=GENO, n_ind=None, n_snp=None):
        self.n_snp = n_snp
        self.mm = np.memmap(path, dtype=np.uint8, mode="r", offset=48,
                            shape=(n_ind, (n_snp + 3) // 4))

    def row_codes(self, i, cols):
        cols = np.asarray(cols)
        byte_idx = cols // 4
        shift = (6 - 2 * (cols % 4)).astype(np.uint8)
        return (np.asarray(self.mm[i, byte_idx]) >> shift) & np.uint8(3)


# --------------------------------------------------------------------------- #
# selection / embedding
# --------------------------------------------------------------------------- #
def select_individuals(ind_ids, anno, lat_win, lon_win):
    """Return (row_idx, meta) for non-modern, dated, geolocated individuals."""
    idx, meta = [], []
    n_modern = n_nocoord = n_nodate = n_window = 0
    for i, gid in enumerate(ind_ids):
        a = anno.get(gid)
        if a is None:
            n_nodate += 1
            continue
        if a["method"].startswith("Modern") or (a["date"] is not None and a["date"] == 0):
            n_modern += 1
            continue
        if a["date"] is None or a["lat"] is None or a["lon"] is None:
            n_nocoord += 1
            continue
        if not (lat_win[0] <= a["lat"] <= lat_win[1]
                and lon_win[0] <= a["lon"] <= lon_win[1]
                and 0 < a["date"] <= 50_000):
            n_window += 1
            continue
        idx.append(i)
        a = dict(a)
        a["_gid"] = gid
        meta.append(a)
    return np.array(idx), meta, dict(modern=n_modern, no_date=n_nodate,
                                     no_coord=n_nocoord, window=n_window)


def pick_sites(chrom, pos, regions, n_cand, seed):
    auto = np.isin(chrom, [str(c) for c in range(1, 23)])
    ok = auto & ~excluded_mask(chrom, pos, regions)
    avail = np.where(ok)[0]
    if len(avail) > n_cand:
        stride = int(np.ceil(len(avail) / n_cand))
        avail = avail[::stride][:n_cand]
    return avail


def maf_screen(geno, ind_idx, cand, n_sub, seed, log=print):
    rng = np.random.default_rng(seed)
    sub = ind_idx if len(ind_idx) <= n_sub else rng.choice(ind_idx, n_sub, replace=False)
    n_pres = np.zeros(len(cand))
    n_alt = np.zeros(len(cand))
    t0 = time.time()
    for i in sub:
        c = geno.row_codes(i, cand)
        n_pres += (c != 3)
        n_alt += (c == 2)
    freq = n_alt / np.maximum(n_pres, 1)
    maf = np.minimum(freq, 1 - freq)
    pres = n_pres / len(sub)
    log(f"  MAF screen over {len(sub)} inds, {len(cand)} cand sites in {time.time()-t0:.1f}s")
    return maf, pres


def build_embedding(geno, site_cols, ind_idx, n_pcs, seed, log=print):
    n_ind, m = len(ind_idx), len(site_cols)
    codes = np.empty((n_ind, m), dtype=np.uint8)
    t0 = time.time()
    for k, i in enumerate(ind_idx):
        codes[k] = geno.row_codes(i, site_cols)
    log(f"  read genotypes {codes.shape} in {time.time()-t0:.1f}s")
    X = np.where(codes == 3, np.float32(np.nan), (codes == 2).astype(np.float32))
    del codes
    col_mean = np.nanmean(X, axis=0)
    ind = np.where(np.isnan(X))
    X[ind] = np.take(col_mean, ind[1])
    X -= X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    X /= sd
    from sklearn.decomposition import PCA
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed)
    T = pca.fit_transform(X).astype(np.float64)
    log(f"  PCA var explained: {np.round(pca.explained_variance_ratio_, 4)}")
    return T


# --------------------------------------------------------------------------- #
# ordered-path trajectories
# --------------------------------------------------------------------------- #
def cluster_bins(T, dates, width, step, seed, kmax, per_cluster, min_cluster):
    """Overlapping time bins -> list of nodes (bin, cluster, centroid, members)."""
    from sklearn.cluster import KMeans

    starts = np.arange(dates.min(), dates.max(), step, dtype=np.float64)
    nodes = []
    for bi, s in enumerate(starts):
        m = np.where((dates >= s) & (dates < s + width))[0]
        if len(m) == 0:
            continue
        if len(m) < min_cluster:
            k = 1
        else:
            k = int(np.clip(len(m) // per_cluster, 1, kmax))
        if k == 1:
            labs = np.zeros(len(m), dtype=int)
        else:
            km = KMeans(n_clusters=k, n_init=3, random_state=seed + bi).fit(T[m])
            labs = km.labels_
        for c in range(k):
            mem = m[labs == c]
            nodes.append({"bin": bi, "start": s, "cluster": c,
                          "centroid": T[mem].mean(axis=0), "members": mem})
    return starts, nodes


def pair_gate(nodes, starts, step, mult):
    """Hard distance gate = mult * median adjacent-bin centroid distance."""
    by_bin = {}
    for n in nodes:
        by_bin.setdefault(n["bin"], []).append(n)
    ds = []
    for b in range(len(starts) - 1):
        A, B = by_bin.get(b, []), by_bin.get(b + 1, [])
        for a in A:
            for c in B:
                ds.append(np.linalg.norm(a["centroid"] - c["centroid"]))
    d0 = float(np.median(ds)) if ds else 0.0
    return mult * d0, d0


def link_chains(nodes, starts, step, gate, max_edges):
    """Greedy one-to-one adjacency linking -> raw chains, split to <= max_edges."""
    by_bin = {}
    for ni, n in enumerate(nodes):
        by_bin.setdefault(n["bin"], []).append(ni)
    edges = []
    for b in range(len(starts) - 1):
        A, B = by_bin.get(b, []), by_bin.get(b + 1, [])
        for a in A:
            for c in B:
                d = np.linalg.norm(nodes[a]["centroid"] - nodes[c]["centroid"])
                if d <= gate:
                    edges.append((d, a, c))
    edges.sort(key=lambda e: e[0])
    succ = {n: -1 for n in range(len(nodes))}
    pred = {n: -1 for n in range(len(nodes))}
    for d, a, c in edges:
        if succ[a] == -1 and pred[c] == -1:
            succ[a] = c
            pred[c] = a
    starts_nodes = [n for n in range(len(nodes)) if pred[n] == -1]
    chains = []
    for s in starts_nodes:
        chain = []
        cur = s
        while cur != -1:
            chain.append(cur)
            cur = succ[cur]
        for j in range(0, len(chain), max_edges + 1):
            chains.append(chain[j:j + max_edges + 1])
    return chains


def keep_segments(chains, nodes, home, width, max_span, min_size):
    """Keep chain segments whose HOME-BIN-assigned member count and span pass."""
    kept = []
    for seg in chains:
        n_assigned = 0
        for nid in seg:
            nb = nodes[nid]["bin"]
            n_assigned += int(sum(1 for m in nodes[nid]["members"] if home[m] == nb))
        span = nodes[seg[-1]]["start"] - nodes[seg[0]]["start"] + width
        if n_assigned >= min_size and span <= max_span:
            kept.append(seg)
    return kept


def home_bins(starts, dates, width):
    centers = starts + width / 2.0
    return np.abs(dates[:, None] - centers[None, :]).argmin(axis=1)


def assign_trajectories(kept, nodes, home, n_ind):
    traj = np.full(n_ind, -1, dtype=int)
    for tid, seg in enumerate(kept):
        for nid in seg:
            n = nodes[nid]
            for m in n["members"]:
                if home[m] == n["bin"]:
                    traj[m] = tid
    return traj


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
def coherence(T, traj):
    from sklearn.metrics import pairwise_distances
    ids = [t for t in np.unique(traj) if t >= 0]
    rng = np.random.default_rng(0)
    within = {}
    for t in ids:
        m = np.where(traj == t)[0]
        if len(m) > 400:
            m = rng.choice(m, 400, replace=False)
        within[t] = float(pairwise_distances(T[m]).mean()) if len(m) > 1 else 0.0
    cross = float(np.median([within[t] for t in ids])) if ids else float("nan")
    return within, cross


def ari(a, b):
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(a, b))


def group_annotation(meta, traj, idx_in_meta):
    """Post-hoc: which published groups fall in which trajectory."""
    out = {}
    for t in sorted(set(traj.tolist())):
        if t < 0:
            continue
        members = idx_in_meta[traj == t]
        groups = [meta[i]["group_id"] for i in members]
        from collections import Counter
        cnt = Counter(groups)
        top = cnt.most_common(5)
        joined = " ".join(g.lower() for g in cnt)
        out[int(t)] = {
            "n": int(len(members)),
            "top_groups": [[g, c] for g, c in top],
            "whg_like": sum(c for g, c in cnt.items() if any(k in g.lower() for k in WHG_KEYS)),
            "eef_like": sum(c for g, c in cnt.items() if any(k in g.lower() for k in EEF_KEYS)),
            "steppe_like": sum(c for g, c in cnt.items() if any(k in g.lower() for k in STEPPE_KEYS)),
        }
    return out


# (kept minimal); trajectory span cap is enforced inside keep_segments.

def run_pipeline(T, dates, width, step, seed, gate_mult, kmax, per_cluster,
                 min_cluster, min_size, max_span, return_nodes=False):
    """Full binning+clustering+linking on the supplied (already embedded) set."""
    starts, nodes = cluster_bins(T, dates, width, step, seed, kmax, per_cluster,
                                 min_cluster)
    gate, d0 = pair_gate(nodes, starts, step, gate_mult)
    max_edges = max(0, int((max_span - width) // step))
    chains = link_chains(nodes, starts, step, gate, max_edges)
    home = home_bins(starts, dates, width)
    kept = keep_segments(chains, nodes, home, width, max_span, min_size)
    traj = assign_trajectories(kept, nodes, home, len(dates))
    if return_nodes:
        return traj, kept, nodes, starts, gate, d0, home
    return traj


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def make_figures(T, dates, traj, meta, idx_in_meta, kept, nodes, metrics, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(outdir, exist_ok=True)
    ids = [t for t in np.unique(traj) if t >= 0]
    cmap = plt.get_cmap("tab20")
    color = {t: cmap(i % 20) for i, t in enumerate(ids)}

    # 1. PCA coloured by trajectory
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    ax[0].scatter(T[:, 0], T[:, 1], s=4, c="#cccccc", label="unassigned")
    for t in ids:
        m = traj == t
        ax[0].scatter(T[m, 0], T[m, 1], s=6, color=color[t], label=f"T{t}")
    ax[0].set(xlabel="PC1", ylabel="PC2",
              title=f"PCA coloured by trajectory ({len(ids)} trajectories)")
    ax[0].legend(fontsize=6, markerscale=2, ncol=2)
    for t in ids:
        m = traj == t
        ax[1].scatter(dates[m], T[m, 0], s=6, color=color[t])
    ax[1].set(xlabel="date (BP)", ylabel="PC1", title="PC1 vs time by trajectory")
    ax[1].invert_xaxis()
    fig.tight_layout()
    fig.savefig(f"{outdir}/traj_v2_pca.png", dpi=130)
    plt.close(fig)

    # 2. date span + spatial extent
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    for i, t in enumerate(ids):
        m = traj == t
        d = dates[m]
        ax[0].plot([d.min(), d.max()], [i, i], "-", color=color[t], lw=3)
        ax[0].scatter([np.median(d)], [i], color="k", s=10, zorder=3)
    ax[0].set(xlabel="date (BP)", ylabel="trajectory (sorted by median date)",
              title="Trajectory date spans (bar) + median (dot)")
    ax[0].invert_xaxis()
    for t in ids:
        m = traj == t
        lon = np.array([meta[j]["lon"] for j in idx_in_meta[m]])
        lat = np.array([meta[j]["lat"] for j in idx_in_meta[m]])
        ax[1].scatter(lon, lat, s=8, color=color[t], alpha=0.6)
    ax[1].set(xlabel="longitude", ylabel="latitude",
              title="Trajectory spatial extent (all trajectories)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/traj_v2_spans.png", dpi=130)
    plt.close(fig)

    # 3. stability
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    boot = metrics["stability"]["bootstrap_ari"]
    ax[0].hist(boot, bins=12, range=(0, 1), color="#4477aa")
    ax[0].axvline(metrics["stability"]["bootstrap_ari_mean"], color="k",
                  label=f"mean={metrics['stability']['bootstrap_ari_mean']:.3f}")
    ax[0].set(xlabel="ARI vs base", ylabel="bootstrap reps",
              title=f"Bootstrap ARI (n={len(boot)})")
    ax[0].legend()
    ks = list(metrics["sensitivity"].keys())
    vs = [metrics["sensitivity"][k] for k in ks]
    ax[1].barh(ks, vs, color="#aa7744")
    ax[1].set(xlim=(0, 1), xlabel="ARI vs base",
              title="Sensitivity to bin width / gate")
    fig.tight_layout()
    fig.savefig(f"{outdir}/traj_v2_stability.png", dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sites", type=int, default=100_000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-pcs", type=int, default=10)
    ap.add_argument("--min-present", type=float, default=0.10,
                    help="coverage floor (fraction of selected sites called)")
    ap.add_argument("--bin-width", type=float, default=1200.0)
    ap.add_argument("--bin-step", type=float, default=400.0)
    ap.add_argument("--gate-mult", type=float, default=1.0)
    ap.add_argument("--kmax", type=int, default=6)
    ap.add_argument("--per-cluster", type=int, default=40)
    ap.add_argument("--min-cluster", type=int, default=10)
    ap.add_argument("--min-traj-size", type=int, default=50)
    ap.add_argument("--max-span", type=float, default=8000.0)
    ap.add_argument("--bootstrap", type=int, default=20)
    ap.add_argument("--maf-frac", type=float, default=0.05)
    ap.add_argument("--pres-frac", type=float, default=0.5)
    ap.add_argument("--out", default="data/derived/trajectories.tsv")
    ap.add_argument("--metrics", default="data/derived/traj_v2_metrics.json")
    ap.add_argument("--fig-prefix", default="fig")
    args = ap.parse_args()

    t_start = time.time()
    log = print
    ind_ids = read_ind()
    anno = read_anno()
    chrom, pos = read_snp()
    log(f"loaded {len(ind_ids)} individuals, {len(pos)} sites")

    ind_idx, meta, drop = select_individuals(ind_ids, anno, (25, 75), (-15, 65))
    log(f"non-modern dated geolocated West Eurasian: {len(ind_idx)} "
        f"(dropped modern={drop['modern']} no_date={drop['no_date']} "
        f"no_coord={drop['no_coord']} out_of_window={drop['window']})")

    n_cand = min(len(pos), max(2 * args.n_sites, 50_000))
    cand = pick_sites(chrom, pos, DEFAULT_EXCLUDE, n_cand, args.seed)
    log(f"candidate autosomal non-excluded sites: {len(cand)}")
    geno = TGENO(n_ind=len(ind_ids), n_snp=len(pos))

    maf, pres = maf_screen(geno, ind_idx, cand, 3000, args.seed, log=log)
    good = (maf >= args.maf_frac) & (maf <= 1 - args.maf_frac) & (pres >= args.pres_frac)
    good_sites = cand[good]
    if len(good_sites) > args.n_sites:
        stride = int(np.ceil(len(good_sites) / args.n_sites))
        good_sites = good_sites[::stride][:args.n_sites]
    log(f"selected sites: {len(good_sites)} (MAF>={args.maf_frac}, pres>={args.pres_frac})")

    # coverage floor + embedding (single genotype pass)
    codes = np.empty((len(ind_idx), len(good_sites)), dtype=np.uint8)
    t0 = time.time()
    for k, i in enumerate(ind_idx):
        codes[k] = geno.row_codes(i, good_sites)
    log(f"  read genotypes {codes.shape} in {time.time()-t0:.1f}s")
    cover = (codes != 3).mean(axis=1)
    keep_cov = cover >= args.min_present
    log(f"coverage floor {args.min_present}: kept {keep_cov.sum()}/{len(keep_cov)}")

    codes = codes[keep_cov]
    meta = [meta[i] for i in np.where(keep_cov)[0]]
    X = np.where(codes == 3, np.float32(np.nan), (codes == 2).astype(np.float32))
    del codes
    col_mean = np.nanmean(X, axis=0)
    ind = np.where(np.isnan(X))
    X[ind] = np.take(col_mean, ind[1])
    X -= X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    X /= sd
    from sklearn.decomposition import PCA
    pca = PCA(n_components=args.n_pcs, svd_solver="randomized", random_state=args.seed)
    T = pca.fit_transform(X).astype(np.float64)
    del X
    log(f"  PCA var explained: {np.round(pca.explained_variance_ratio_, 4)}")
    dates = np.array([a["date"] for a in meta], dtype=np.float64)
    n = len(dates)

    traj, kept, nodes, starts, gate, d0, home = run_pipeline(
        T, dates, args.bin_width, args.bin_step, args.seed, args.gate_mult,
        args.kmax, args.per_cluster, args.min_cluster, args.min_traj_size,
        args.max_span, return_nodes=True)
    ids = [t for t in np.unique(traj) if t >= 0]
    log(f"trajectories kept: {len(ids)} over {n} individuals; "
        f"gate={gate:.4f} (median adjacent centroid dist d0={d0:.4f}); "
        f"{int((traj < 0).sum())} unassigned")

    # ---- per-trajectory stats + coherence ----
    within, cross = coherence(T, traj)
    idx_in_meta = np.arange(n)
    per_traj = {}
    max_span_seen, min_size_seen = 0.0, 10 ** 9
    for t in ids:
        m = traj == t
        d = dates[m]
        lons = np.array([meta[j]["lon"] for j in np.where(m)[0]])
        lats = np.array([meta[j]["lat"] for j in np.where(m)[0]])
        span = float(d.max() - d.min())
        size = int(m.sum())
        max_span_seen = max(max_span_seen, span)
        min_size_seen = min(min_size_seen, size)
        per_traj[int(t)] = {
            "n": size, "date_min": float(d.min()), "date_max": float(d.max()),
            "date_median": float(np.median(d)), "span_years": span,
            "lat_min": float(lats.min()), "lat_max": float(lats.max()),
            "lon_min": float(lons.min()), "lon_max": float(lons.max()),
            "coherence_mean_pairwise": within[int(t)],
            "coherence_ratio_vs_cross": within[int(t)] / cross if cross else None,
        }
    # assert caps
    assert max_span_seen <= args.max_span + 1e-6, \
        f"span cap violated: {max_span_seen} > {args.max_span}"
    assert min_size_seen >= args.min_traj_size, \
        f"min size violated: {min_size_seen} < {args.min_traj_size}"
    log(f"CHECK span cap: max observed span={max_span_seen:.0f} <= {args.max_span:.0f} OK")
    log(f"CHECK min size: min observed size={min_size_seen} >= {args.min_traj_size} OK")

    gro = group_annotation(meta, traj, idx_in_meta)

    # ---- bootstrap ----
    rng = np.random.default_rng(args.seed + 1000)
    aris = []
    t0 = time.time()
    for b in range(args.bootstrap):
        pick = np.unique(rng.integers(0, n, size=n))
        tb = run_pipeline(T[pick], dates[pick], args.bin_width, args.bin_step,
                          args.seed + 1 + b, args.gate_mult, args.kmax,
                          args.per_cluster, args.min_cluster, args.min_traj_size,
                          args.max_span)
        aris.append(ari(traj[pick], tb))
    bar = float(np.mean(aris)); bsd = float(np.std(aris))
    log(f"BOOTSTRAP ({len(aris)} reps) ARI mean={bar:.3f} sd={bsd:.3f} "
        f"({time.time()-t0:.1f}s)")

    # ---- sensitivity ----
    sens = {}
    t_w = run_pipeline(T, dates, args.bin_width + args.bin_step, args.bin_step,
                       args.seed, args.gate_mult, args.kmax, args.per_cluster,
                       args.min_cluster, args.min_traj_size, args.max_span)
    sens["bin_width+400"] = ari(traj, t_w)
    t_g = run_pipeline(T, dates, args.bin_width, args.bin_step, args.seed,
                       args.gate_mult * 1.5, args.kmax, args.per_cluster,
                       args.min_cluster, args.min_traj_size, args.max_span)
    sens["gate_x1.5"] = ari(traj, t_g)
    t_g2 = run_pipeline(T, dates, args.bin_width, args.bin_step, args.seed,
                        args.gate_mult * 0.7, args.kmax, args.per_cluster,
                        args.min_cluster, args.min_traj_size, args.max_span)
    sens["gate_x0.7"] = ari(traj, t_g2)
    t_s = run_pipeline(T, dates, args.bin_width, args.bin_step * 2, args.seed,
                       args.gate_mult, args.kmax, args.per_cluster,
                       args.min_cluster, args.min_traj_size, args.max_span)
    sens["bin_step800"] = ari(traj, t_s)
    log(f"SENSITIVITY ARI: { {k: round(v,3) for k,v in sens.items()} }")

    # ---- write table ----
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    hdr = ["individual", "date_BP", "lat", "lon", "group_id"] + \
          [f"PC{i+1}" for i in range(T.shape[1])] + ["trajectory_id", "bin_id"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for j in range(n):
            row = [meta[j]["_gid"], f"{dates[j]:.0f}",
                   f"{meta[j]['lat']:.5f}", f"{meta[j]['lon']:.5f}",
                   meta[j]["group_id"]] + \
                  [f"{T[j, c]:.5f}" for c in range(T.shape[1])] + \
                  [str(traj[j]), str(int(home[j]))]
            fh.write("\t".join(row) + "\n")
    log(f"wrote {args.out}")

    metrics = {
        "config": vars(args),
        "counts": {"n_ind_total": len(ind_ids), "n_sites_total": len(pos),
                   "dropped": drop, "n_usable": int(len(ind_idx)),
                   "n_sites_candidate": int(len(cand)),
                   "n_sites_selected": int(len(good_sites)),
                   "n_coverage_kept": int(keep_cov.sum()),
                   "n_bins": int(len(starts)), "n_nodes": int(len(nodes)),
                   "n_trajectories": len(ids),
                   "n_unassigned": int((traj < 0).sum())},
        "pca_var_explained": [float(v) for v in pca.explained_variance_ratio_],
        "gate": {"value": float(gate), "median_adjacent_centroid_dist": float(d0),
                 "gate_mult": args.gate_mult},
        "span_cap_check": {"max_span_seen": max_span_seen, "cap": args.max_span},
        "min_size_check": {"min_size_seen": min_size_seen, "floor": args.min_traj_size},
        "trajectories": per_traj,
        "cross_trajectory_mean_pairwise": cross,
        "group_annotation": gro,
        "stability": {"bootstrap_ari": [float(a) for a in aris],
                      "bootstrap_ari_mean": bar, "bootstrap_ari_sd": bsd,
                      "n_bootstrap": len(aris)},
        "sensitivity": {k: float(v) for k, v in sens.items()},
    }
    # add genetic id into meta for the table (was not stored)
    with open(args.metrics, "w") as fh:
        json.dump(metrics, fh, indent=2)
    log(f"wrote {args.metrics}")

    make_figures(T, dates, traj, meta, idx_in_meta, kept, nodes, metrics,
                 args.fig_prefix)
    log(f"wrote figures traj_v2_*.png")
    log(f"TOTAL RUNTIME {time.time()-t_start:.1f}s")


if __name__ == "__main__":
    main()
