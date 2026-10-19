#!/usr/bin/env python3
"""omega-style statistic, fixed: whole-branch length vs a genome-wide leaf-count
expectation, plus a generalised-skyline Ne(t) for the demographic part.

Why the first attempt failed (docs/arg-metric-ground-rules): the expectation was a mean
RAW length per (band, leaf count) taken locally, and the score subtracted "bands touched",
so a branch entering a wide band briefly scored ~0 (LCT gave omega = -1.0).

Fixed construction:
  * per branch b: total length L_b = t[parent] - t[node]   (whole branch, no per-band split)
  * leaf count n_b = number of descendant samples
  * expectation E(n) = median of L over ALL branches genome-wide (sampled) with the same
    log2 leaf-count bin   -> demography-free, like the clade-size-matched null
  * omega_b = L_b / E(n_b);  per-SNP score = omega_b - 1
    (>0 means the branch is longer than branches carrying as many samples, i.e. a long
    branch above a clade that coalesced recently = the sweep signature)

Modes
  null <chrom> [stride]        build E(n) from sampled trees of a chromosome -> npz
  scan <chrom> <lo> <hi>       per-SNP score distribution in a window + named loci
  skyline <chrom> [stride]     generalised skyline Ne(t) from the ARG
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
POPFILE = OUT / 'samples_population.txt'
BANDS = np.array([0, 20, 100, 500, 2000, 10_000, 50_000, 200_000, 1e7], dtype=float)


def eur_haps():
    sel = []
    with open(POPFILE) as f:
        next(f)
        for i, line in enumerate(f):
            p = line.split()
            if len(p) >= 7 and p[6] == 'EUR':
                sel.append(i)
    return [h for i in sel for h in (2 * i, 2 * i + 1)]


def chunks_for(chrom):
    return sorted(OUT.glob(f'chr{chrom}_*.dated.trees'))


def leaf_counts(tree, t, n_nodes, samples):
    cnt = np.zeros(n_nodes, dtype=np.int32)
    cnt[samples] = 1
    par = tree.parent_array
    for v in np.argsort(t, kind='stable'):
        p = par[v]
        if p != -1:
            cnt[p] += cnt[v]
    return cnt


def branches_in(ts, stride=1, pos_range=None, start=0):
    """Yield (leafbin, length, abs_pos_of_tree_mid, child) for branches in the range."""
    t = ts.tables.nodes.time
    samples = ts.samples()
    out = []
    for ti, tree in enumerate(ts.trees()):
        if stride > 1 and ti % stride:
            continue
        mid = tree.interval.left + start
        if pos_range is not None and not (pos_range[0] <= mid < pos_range[1]):
            continue
        cnt = leaf_counts(tree, t, ts.num_nodes, samples)
        par = tree.parent_array
        kids = np.where(par != -1)[0]
        for c in kids:
            # same zero-length collapse as scan(), so null and observations are commensurate
            v = int(c); pp = int(par[v])
            while pp != -1 and (t[pp] - t[v]) < 1e-3:
                v = pp; pp = int(par[v])
            if pp == -1 or not np.isfinite(t[pp]):
                continue
            out.append((int(np.floor(np.log2(max(int(cnt[c]), 1)))), float(t[pp] - t[c]), float(mid), int(c)))
    return out


def build_null(chrom, stride=100):
    rows = []
    for f in chunks_for(chrom):
        ts = tskit.load(str(f)).simplify(samples=eur_haps())
        start = int(f.name.split('_')[1])
        rows.extend(branches_in(ts, stride=stride, start=start))
    df = pd.DataFrame(rows, columns=['leafbin', 'length', 'mid', 'child'])
    E = df.groupby('leafbin').length.agg(['median', 'mean', 'size'])
    E.to_csv(OUT / f'omega_null_chr{chrom}.csv')
    print(f'chr{chrom}: null from {len(df):,} sampled branches')
    print(E.to_string())


def scan(chrom, lo, hi, targets):
    """Per-SNP omega in a window, using the chromosome null (built genome-wide from chr files)."""
    e = {}
    for f in sorted(OUT.glob('omega_null_chr*.csv')):
        d = pd.read_csv(f)
        for r in d.itertuples(index=False):
            e.setdefault(int(r[0]), []).append(float(r[1]))        # leafbin, median
    E = {k: float(np.median(v)) for k, v in e.items()}
    rows = []
    for f in chunks_for(chrom):
        start = int(f.name.split('_')[1]); end = int(f.name.split('_')[2].split('.')[0])
        if end <= lo or start >= hi:
            continue
        ts = tskit.load(str(f)).simplify(samples=eur_haps())
        t = ts.tables.nodes.time
        samples = ts.samples()
        for tree in ts.trees():
            l_, r_ = tree.interval.left + start, tree.interval.right + start
            if r_ <= lo or l_ >= hi:
                continue
            cnt = leaf_counts(tree, t, ts.num_nodes, samples)
            par = tree.parent_array
            Lb = {}
            for c in np.where(par != -1)[0]:
                # collapse zero-length chains: the allele's lineage continues upward
                v = int(c); pp = int(par[v])
                while pp != -1 and (t[pp] - t[v]) < 1e-3:   # collapse degenerate branches
                    v = pp; pp = int(par[v])
                if pp == -1 or not np.isfinite(t[pp]):
                    continue
                Lb[int(c)] = (float(t[pp] - t[c]),
                              int(np.floor(np.log2(max(int(cnt[c]), 1)))))
            for m in tree.mutations():
                pp = ts.site(m.site).position + start
                if not (lo <= pp < hi):
                    continue
                if int(m.node) not in Lb:
                    continue
                L, lb = Lb[int(m.node)]
                exp = E.get(lb)
                if not exp:
                    continue
                rows.append(dict(pos=int(pp), leafbin=lb, length=L, omega=L / exp))
    df = pd.DataFrame(rows)
    if not len(df):
        print('no SNPs'); return
    df['score'] = df.omega - 1
    print(f'chr{chrom}:{lo}-{hi} | SNPs {len(df):,} | median omega {df.omega.median():.3f} '
          f'| median score {df.score.median():.3f}')
    for lab, p in targets:
        m = df[np.abs(df.pos - p) <= 1]
        if not len(m):
            m = df[np.abs(df.pos - p) <= 1000]
        if not len(m):
            print(f'  {lab}: absent'); continue
        r = m.iloc[np.abs(m.pos - p).argmin()]
        pct = 100 * (df.score < r.score).mean()
        print(f'  {lab:10s} pos={int(r.pos)} leafbin={int(r.leafbin)} length={r.length:.0f} gen '
              f'omega={r.omega:.2f} score={r.score:+.2f} percentile={pct:.1f}%')


def skyline(chrom, stride=100):
    """Generalised skyline: Ne_j = k_j(k_j-1)dt_j / (4 c_j)."""
    k_acc = [[] for _ in BANDS[:-1]]
    c_acc = [0 for _ in BANDS[:-1]]
    for f in chunks_for(chrom):
        ts = tskit.load(str(f)).simplify(samples=eur_haps())
        t = ts.tables.nodes.time
        for ti, tree in enumerate(ts.trees()):
            if stride > 1 and ti % stride:
                continue
            par = tree.parent_array
            kids = np.where(par != -1)[0]
            for bi in range(len(BANDS) - 1):
                a, b = BANDS[bi], BANDS[bi + 1]
                lin = 0
                for c in kids:
                    tc, tp = t[c], t[par[c]]
                    if tc < b and tp > a:
                        lin += 1
                if lin:
                    k_acc[bi].append(lin)
            # coalescent events = internal nodes
            for v in np.where(t > 0)[0]:
                if par[v] == -1:
                    continue
                tv = t[v]
                for bi in range(len(BANDS) - 1):
                    if BANDS[bi] <= tv < BANDS[bi + 1]:
                        c_acc[bi] += 1
                        break
    print('band\twidth\tmean_lineages\tcoalescences\tNe_est')
    for bi in range(len(BANDS) - 1):
        dt = BANDS[bi + 1] - BANDS[bi]
        if not k_acc[bi] or c_acc[bi] == 0:
            continue
        k = float(np.mean(k_acc[bi]))
        ne = k * (k - 1) * dt / (4 * c_acc[bi])
        print(f'{BANDS[bi]:.0f}-{BANDS[bi+1]:.0f}\t{dt:.0f}\t{k:.1f}\t{c_acc[bi]}\t{ne:.0f}')


def main():
    mode = sys.argv[1]
    chrom = int(sys.argv[2])
    if mode == 'null':
        build_null(chrom, int(sys.argv[3]) if len(sys.argv) > 3 else 100)
    elif mode == 'skyline':
        skyline(chrom, int(sys.argv[3]) if len(sys.argv) > 3 else 100)
    elif mode == 'scan':
        lo, hi = int(sys.argv[3]), int(sys.argv[4])
        targets = [(x.split(':')[0], int(x.split(':')[1])) for x in sys.argv[5].split(',')] if len(sys.argv) > 5 else []
        scan(chrom, lo, hi, targets)


if __name__ == '__main__':
    main()
