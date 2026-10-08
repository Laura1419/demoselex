#!/usr/bin/env python3
"""Stage A prototype: infer temporal population trajectories from genome-wide
genetic similarity plus time (docs/D-method.md, Section 2, primary route).

Method
------
1. Read a bounded, MAF-filtered, low-missingness subset of autosomal sites from
   the AADR 1240K EIGENSTRAT triplet.  A list of *excluded* chromosome regions
   is always applied, so any candidate locus (e.g. chr6 MHC, chr1 AMY1/AMY2)
   can be removed from the embedding before trajectory inference.
2. Build a low-dimensional PCA embedding on genotypes (mean-imputed,
   variance-standardised), random seed fixed.
3. Cluster with a *time-aware* constraint: samples may be linked only if their
   radiocarbon dates differ by at most ``--temporal-window`` years, OR if the
   pair are reciprocal genetic nearest neighbours.  A symmetric genetic kNN
   graph over the PCA embedding is pruned by that rule, then partitioned into
   ``--k`` trajectories with spectral clustering.  Trajectories may start/end
   (disconnected in the gated graph) and are not named after archaeology.
4. Stability: the run is repeated under a second parameter setting and a
   per-individual bootstrap; agreement is reported as an adjusted Rand index.

The AADR v66 ``.geno`` file used here is the transposed "TGENO" layout:
48-byte header, then ``n_ind`` consecutive blocks of ``ceil(n_snp/4)`` bytes,
each byte holding four consecutive SNPs' 2-bit genotype codes (0,2 = the two
homozygous states, 3 = missing; code 1 is absent in this patched release).
"""
from __future__ import annotations

import argparse
import time

import numpy as np

GENO = "data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.geno"
IND = "data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.ind"
SNP = "data/raw/aadr-v66.p1/v66.p1_1240K.aadr.patch.PUB.snp"
ANNO = "data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.anno"

# GRCh37 coordinates.  Default exclusion: chr6 MHC and chr1 AMY1/AMY2.
DEFAULT_EXCLUDE = (("6", 28_500_000, 33_500_000), ("1", 103_500_000, 104_500_000))


# --------------------------------------------------------------------------- #
# metadata / site selection
# --------------------------------------------------------------------------- #
def read_ind(path=IND):
    ids = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if parts:
                ids.append(parts[0])
    return ids


def read_anno(path=ANNO):
    """Genetic ID -> (date_BP, lat, lon).  First row per ID wins."""
    import csv

    info = {}
    with open(path, newline="") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr)
        for row in rdr:
            gid = row[0]
            if gid in info:
                continue
            try:
                date = float(row[10])
                lat = float(row[17])
                lon = float(row[18])
            except (ValueError, IndexError):
                continue
            info[gid] = (date, lat, lon)
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
    """Boolean mask: True where a site falls inside an excluded region."""
    bad = np.zeros(len(pos), dtype=bool)
    for c, lo, hi in regions:
        bad |= (chrom == c) & (pos >= lo) & (pos < hi)
    return bad


def candidate_sites(chrom, pos, regions, n_cand, seed):
    """Autosomal, non-excluded sites evenly spaced across the genome."""
    auto = np.isin(chrom, [str(c) for c in range(1, 23)])
    ok = auto & ~excluded_mask(chrom, pos, regions)
    avail = np.where(ok)[0]
    if len(avail) > n_cand:
        stride = int(np.ceil(len(avail) / n_cand))
        avail = avail[::stride][:n_cand]
    return avail


# --------------------------------------------------------------------------- #
# genotype I/O (transposed TGENO layout)
# --------------------------------------------------------------------------- #
class TGENO:
    def __init__(self, path=GENO, n_ind=None, n_snp=None):
        self.n_ind = n_ind
        self.n_snp = n_snp
        self.bpi = (n_snp + 3) // 4  # bytes per individual
        self.mm = np.memmap(path, dtype=np.uint8, mode="r", offset=48,
                            shape=(n_ind, self.bpi))

    def row_codes(self, i, cols=None):
        """2-bit codes for individual ``i`` at site columns ``cols``.

        Returns uint8 array; 3 marks missing.  With ``cols=None`` all SNPs.
        """
        if cols is None:
            b = np.asarray(self.mm[i])
            c = np.empty(b.size * 4, dtype=np.uint8)
            c[0::4] = (b >> 6) & 3
            c[1::4] = (b >> 4) & 3
            c[2::4] = (b >> 2) & 3
            c[3::4] = b & 3
            return c[: self.n_snp]
        cols = np.asarray(cols)
        byte_idx = cols // 4
        shift = (6 - 2 * (cols % 4)).astype(np.uint8)
        b = np.asarray(self.mm[i, byte_idx])
        return (b >> shift) & 3


# --------------------------------------------------------------------------- #
# embedding
# --------------------------------------------------------------------------- #
def build_embedding(geno, site_cols, ind_idx, n_pcs, seed, log=print):
    """PCA embedding (ind_idx x n_pcs) over ``site_cols``.

    Missing genotypes are mean-imputed, then each site is variance-standardised.
    """
    n_ind = len(ind_idx)
    m = len(site_cols)
    X = np.empty((n_ind, m), dtype=np.float32)
    t0 = time.time()
    for k, i in enumerate(ind_idx):
        codes = geno.row_codes(i, cols=site_cols)
        d = np.where(codes == 3, np.nan, (codes == 2).astype(np.float32))
        X[k] = d
    log(f"  read genotypes {X.shape} in {time.time()-t0:.1f}s")

    col_mean = np.nanmean(X, axis=0)
    ind = np.where(np.isnan(X))
    X[ind] = np.take(col_mean, ind[1])
    col_sd = X.std(axis=0)
    col_sd[col_sd == 0] = 1.0
    X /= col_sd
    X -= X.mean(axis=0)

    from sklearn.decomposition import PCA

    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed)
    T = pca.fit_transform(X)
    log(f"  PCA var explained: {np.round(pca.explained_variance_ratio_, 3)}")
    return T.astype(np.float64)


# --------------------------------------------------------------------------- #
# time-aware clustering
# --------------------------------------------------------------------------- #
def time_aware_affinity(T, dates, window, gnn_big=30, gnn_small=15):
    """Symmetric normalised affinity with a temporal gate.

    An edge survives if |date_i - date_j| <= window, or the two samples are
    reciprocal members of each other's ``gnn_small`` genetic nearest neighbours.
    """
    from sklearn.neighbors import NearestNeighbors
    from scipy.sparse import csr_matrix, coo_matrix

    n = len(dates)
    nn = NearestNeighbors(n_neighbors=min(gnn_big, n), metric="euclidean").fit(T)
    dist, idx = nn.kneighbors(T)

    rows = np.repeat(np.arange(n), gnn_big - 1)
    cols = idx[:, 1:].ravel()
    up = rows < cols
    a, b = rows[up], cols[up]

    # reciprocal small-neighbour test (vectorised via encoded pair ids)
    gs = min(gnn_small, n) - 1
    sn = idx[:, 1:1 + gs]
    codes_small = (np.repeat(np.arange(n), gs).astype(np.int64) * n
                   + sn.ravel())
    Ss = np.unique(codes_small)
    codes_ab = a.astype(np.int64) * n + b
    codes_ba = b.astype(np.int64) * n + a
    recip = np.isin(codes_ab, Ss) & np.isin(codes_ba, Ss)

    dt = np.abs(dates[a] - dates[b])
    keep = (dt <= window) | recip
    a, b = a[keep], b[keep]

    sigma = float(np.median(dist[:, 1])) + 1e-9
    d = np.linalg.norm(T[a] - T[b], axis=1)
    w = np.exp(-(d * d) / (2 * sigma * sigma))
    w = np.maximum(w, 1e-4)
    rows2 = np.concatenate([a, b])
    cols2 = np.concatenate([b, a])
    vals = np.concatenate([w, w])
    W = coo_matrix((vals, (rows2, cols2)), shape=(n, n)).tocsr()
    W = W.maximum(W.T)
    deg = np.asarray(W.sum(axis=1)).ravel()
    deg[deg == 0] = 1.0
    Dm = csr_matrix((1.0 / np.sqrt(deg), (np.arange(n), np.arange(n))),
                    shape=(n, n))
    return (Dm @ W @ Dm).tocsr()


def cluster_time_aware(T, dates, k, window, seed):
    from sklearn.cluster import SpectralClustering

    W = time_aware_affinity(T, dates, window)
    sc = SpectralClustering(n_clusters=k, affinity="precomputed",
                            assign_labels="kmeans", random_state=seed)
    labels = sc.fit_predict(W)
    return labels


def ari(a, b):
    from sklearn.metrics import adjusted_rand_score

    return adjusted_rand_score(a, b)


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def individual_filter(geno, site_cols, ind_idx, min_frac, log=print):
    """Keep individuals with at least ``min_frac`` of sites genotyped."""
    t0 = time.time()
    present = np.empty(len(ind_idx))
    for k, i in enumerate(ind_idx):
        codes = geno.row_codes(i, cols=site_cols)
        present[k] = (codes != 3).mean()
    mask = present >= min_frac
    log(f"  individual coverage filter: {mask.sum()}/{len(mask)} kept "
        f"(min_frac={min_frac}) in {time.time()-t0:.1f}s")
    return mask


def load_matrix(n_sites, seed, n_cand_ind=4000, min_present=0.3, log=print):
    ind_ids = read_ind()
    anno = read_anno()
    chrom, pos = read_snp()

    keep_meta = []
    for i, gid in enumerate(ind_ids):
        if gid in anno:
            date, lat, lon = anno[gid]
            if 25 <= lat <= 75 and -15 <= lon <= 65:
                keep_meta.append((i, gid, date))
    log(f"West Eurasian individuals with dates: {len(keep_meta)}")

    n_cand = min(len(pos), max(4 * n_sites, 50_000))
    cand = candidate_sites(chrom, pos, DEFAULT_EXCLUDE, n_cand, seed)
    log(f"candidate sites after region exclusion: {len(cand)}")

    geno = TGENO(n_ind=len(ind_ids), n_snp=len(pos))

    # pass 1: MAF / missingness from a random subsample of individuals
    rng = np.random.default_rng(seed)
    sub = np.array([m[0] for m in keep_meta])
    if len(sub) > n_cand_ind:
        sub = rng.choice(sub, size=n_cand_ind, replace=False)
    n_pres = np.zeros(len(cand))
    n_alt = np.zeros(len(cand))
    t0 = time.time()
    for i in sub:
        codes = geno.row_codes(i, cols=cand)
        n_pres += (codes != 3)
        n_alt += (codes == 2)
    freq = n_alt / np.maximum(n_pres, 1)
    maf = np.minimum(freq, 1 - freq)
    pres_frac = n_pres / len(sub)
    log(f"  pass1 MAF over {len(sub)} inds in {time.time()-t0:.1f}s")

    good = (maf >= 0.05) & (maf <= 0.95) & (pres_frac >= 0.5)
    good_cand = cand[good]
    if len(good_cand) > n_sites:
        stride = int(np.ceil(len(good_cand) / n_sites))
        good_cand = good_cand[::stride][:n_sites]
    log(f"selected sites: {len(good_cand)} (MAF>=0.05, pres>=0.5)")

    # pass 2: build embedding matrix for all kept individuals
    ind_idx = np.array([m[0] for m in keep_meta])
    dates = np.array([m[2] for m in keep_meta], dtype=np.float64)
    fmask = individual_filter(geno, good_cand, ind_idx, min_present, log=log)
    ind_idx, dates = ind_idx[fmask], dates[fmask]
    T = build_embedding(geno, good_cand, ind_idx, 10, seed, log=log)
    return T, dates, ind_ids, ind_idx, good_cand, maf[good][:len(good_cand)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sites", type=int, default=100_000)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--temporal-window", type=float, default=2000.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-pcs", type=int, default=10)
    ap.add_argument("--out", default="data/derived/trajectories.tsv")
    ap.add_argument("--fig", default="fig/traj_pca_time.png")
    ap.add_argument("--n-inds", type=int, default=4000,
                    help="subsample used for MAF estimation")
    ap.add_argument("--bootstrap", type=int, default=5)
    ap.add_argument("--min-present", type=float, default=0.3)
    ap.add_argument("--n-inds-maf", type=int, default=3000)
    args = ap.parse_args()

    t_start = time.time()
    log = print
    T, dates, ind_ids, ind_idx, site_cols, maf_sel = load_matrix(
        args.n_sites, args.seed, n_cand_ind=args.n_inds,
        min_present=args.min_present, log=log)
    n = len(dates)
    labels = cluster_time_aware(T, dates, args.k, args.temporal_window, args.seed)
    log(f"k={args.k} window={args.temporal_window:.0f} seed={args.seed}: "
        f"{len(np.unique(labels))} trajectories over {n} individuals")

    # ---- stability: alternative parameters ----
    k2 = args.k + 2
    w2 = max(500.0, args.temporal_window / 2)
    seed2 = args.seed + 1
    labels2 = cluster_time_aware(T, dates, k2, w2, seed2)
    ari_params = ari(labels, labels2)
    log(f"STABILITY param ARI (k={args.k},w={args.temporal_window:.0f},s={args.seed}) "
        f"vs (k={k2},w={w2:.0f},s={seed2}) = {ari_params:.3f}")

    # ---- bootstrap over individuals ----
    rng = np.random.default_rng(args.seed + 100)
    aris = []
    for b in range(args.bootstrap):
        pick = rng.integers(0, n, size=n)
        uniq = np.unique(pick)
        Tb = T[uniq]
        try:
            lb = cluster_time_aware(Tb, dates[uniq], args.k, args.temporal_window,
                                    args.seed + b)
            aris.append(ari(labels[uniq], lb))
        except Exception as exc:  # pragma: no cover - degenerate resample
            log(f"  bootstrap {b} skipped: {exc}")
    boot_mean = float(np.mean(aris)) if aris else float("nan")
    boot_sd = float(np.std(aris)) if aris else float("nan")
    log(f"BOOTSTRAP ({len(aris)} reps) mean ARI={boot_mean:.3f} sd={boot_sd:.3f}")

    # ---- write table ----
    import os
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    hdr = ["individual", "date_BP"] + [f"PC{i+1}" for i in range(T.shape[1])] + \
          ["trajectory_id", "k", "temporal_window", "seed", "n_sites"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for j in range(n):
            row = [ind_ids[ind_idx[j]], f"{dates[j]:.0f}"] + \
                  [f"{T[j, c]:.5f}" for c in range(T.shape[1])] + \
                  [str(labels[j]), str(args.k), f"{args.temporal_window:.0f}",
                   str(args.seed), str(len(site_cols))]
            fh.write("\t".join(row) + "\n")
    log(f"wrote {args.out}")

    # ---- figure ----
    make_figure(T, dates, labels, args.fig, args.k, args.temporal_window)
    log(f"wrote {args.fig}")
    log(f"total runtime {time.time()-t_start:.1f}s")


def make_figure(T, dates, labels, path, k, window):
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    cmap = plt.get_cmap("tab10")
    for lab in np.unique(labels):
        m = labels == lab
        axes[0].scatter(T[m, 0], T[m, 1], s=6, color=cmap(lab % 10),
                        label=f"traj {lab}", alpha=0.7)
        axes[1].scatter(dates[m], T[m, 0], s=6, color=cmap(lab % 10), alpha=0.7)
        axes[2].scatter(T[m, 1], T[m, 2], s=6, color=cmap(lab % 10), alpha=0.7)
    axes[0].set(xlabel="PC1", ylabel="PC2",
                title=f"PC1-PC2, {k} trajectories (window={window:.0f}y)")
    axes[0].legend(fontsize=7, markerscale=2)
    axes[1].set(xlabel="date (BP)", ylabel="PC1", title="PC1 vs time")
    axes[2].set(xlabel="PC2", ylabel="PC3", title="PC2-PC3")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
