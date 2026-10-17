#!/usr/bin/env python3
"""Composite null: clade size x local recombination, pooled genome-wide.

Motivation (from three failed nulls):
  local-scope null      -> absorbs the sweep itself
  size-only null        -> regionally miscalibrated (telomeres dominate, LCT invisible)
  size x local-rho null -> controls local demography via rho, keeps the sweep visible

rho for a position = number of tree breakpoints within +/-250 kb, per bp.
Clades are binned by log2(descendant sample count) x rho bin; the null for a variant is
the age distribution of clades in the same (size bin +/- 1, rho bin +/- 1) cells.

Stages:
  build-null <every_nth>   sample chunks genome-wide, collect clade (size,rho,age) cells
  score <chrom>            per-variant p_age for the chromosome
  region                   regional excess test (1 Mb windows) + BH FDR + telomere mask

Writes data/derived/arggenome/composite_*.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/argphased'
HALF = 250_000           # rho window half-width (bp)
SEGDUP = ROOT / 'data/raw/annotation/segdup/genomicSuperDups.txt.gz'
_MASK = None


def segdup_mask():
    """Merged segdup intervals per chromosome: {chrom: (starts, ends)}."""
    global _MASK
    if _MASK is not None:
        return _MASK
    import gzip
    from collections import defaultdict
    iv = defaultdict(list)
    with gzip.open(SEGDUP, 'rt') as f:
        for line in f:
            p = line.split('\t')
            if len(p) < 4 or not p[1].startswith('chr'):
                continue
            c = p[1].replace('chr', '')
            if not c.isdigit() or not 1 <= int(c) <= 22:
                continue
            iv[int(c)].append((int(p[2]), int(p[3])))
    out = {}
    for c, lst in iv.items():
        lst.sort()
        merged = []
        for a, b in lst:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        out[c] = (np.array([m[0] for m in merged]), np.array([m[1] for m in merged]))
    _MASK = out
    print(f'segdup mask: {sum(len(v[0]) for v in out.values())} merged intervals on {len(out)} chromosomes',
          flush=True)
    return out


def is_segdup(chrom, pos):
    m = segdup_mask().get(chrom)
    if m is None:
        return False
    starts, ends = m
    i = np.searchsorted(starts, pos, side='right') - 1
    return bool(i >= 0 and pos <= ends[i])
RHO_EDGES = np.array([0.0, 0.05, 0.1, 0.2, 0.4, 0.8, np.inf])   # breakpoints per kb


def rho_at(bounds, pos, half=HALF):
    """Breakpoints per kb within +/-half bp of pos (bounds = sorted breakpoint positions)."""
    lo = np.searchsorted(bounds, pos - half)
    hi = np.searchsorted(bounds, pos + half)
    return (hi - lo) / (2 * half / 1000.0)


def rbin_of(r):
    return int(np.searchsorted(RHO_EDGES, r, side='right') - 1)


def chunk_files(chrom=None):
    pat = f'chr{chrom}_*_*.dated.trees' if chrom else 'chr*_*_*.dated.trees'
    return sorted(OUT.glob(pat))


def build_null(every=5, chrom=None, part=None, tree_stride=5):
    files = chunk_files(chrom)[::every]
    cells = {}
    for f in files:
        chrom = int(f.name.split('_')[0].replace('chr', ''))
        start = int(f.name.split('_')[1])
        ts = tskit.load(str(f))
        t = ts.tables.nodes.time
        samples = ts.samples()
        bounds = np.array(sorted({int(f.name.split('_')[1])} |
                                 {int(iv.left) + start for iv in ts.breakpoints()})) if False else None
        # tree breakpoints in absolute coords
        bps = np.array(list(ts.breakpoints())) + start
        for ti, tree in enumerate(ts.trees()):
            if tree_stride > 1 and ti % tree_stride:
                continue
            par = tree.parent_array
            counts = np.zeros(ts.num_nodes, dtype=np.int32)
            for s in samples:
                v = s
                while v != -1:
                    counts[v] += 1
                    v = par[v]
            nz = np.nonzero(counts)[0]
            if len(nz) == 0:
                continue
            mid = (tree.interval.left + tree.interval.right) / 2 + start
            rb = rbin_of(rho_at(bps, mid))
            comp = 1 if is_segdup(chrom, int(mid)) else 0
            sb = np.round(np.log2(np.clip(counts[nz], 1, None))).astype(int)
            for s_bin in np.unique(sb):
                key = (int(s_bin), rb, comp)
                cells.setdefault(key, []).append(t[nz[sb == s_bin]])
        print(f'  {f.name}: cells so far {len(cells)}', flush=True)
    arr = {f'{k[0]}_{k[1]}_{k[2]}': np.sort(np.concatenate(v)) for k, v in cells.items()}
    name = f'composite_null_{part}.npz' if part else 'composite_null.npz'
    np.savez_compressed(OUT / name, **arr)
    print(f'wrote {name} with {len(arr)} cells, '
          f'{sum(len(v) for v in arr.values()):,} clades')


def score(chrom):
    z = np.load(OUT / 'composite_null.npz')
    rows = []
    for f in chunk_files(chrom):
        start = int(f.name.split('_')[1])
        ts = tskit.load(str(f))
        t = ts.tables.nodes.time
        G = ts.genotype_matrix()
        spos = np.array([s.position for s in ts.sites()]) + start
        bps = np.array(list(ts.breakpoints())) + start
        mutnode = {}
        for m in ts.mutations():
            p_abs = int(ts.site(m.site).position) + start
            mutnode.setdefault(p_abs, m.node)
        for j, p_abs in enumerate(spos):
            carriers = int(G[j].sum())
            if carriers < 3:
                continue
            node = mutnode.get(int(p_abs))
            if node is None:
                continue
            age = float(t[node])
            sb = int(round(np.log2(carriers)))
            rb = rbin_of(rho_at(bps, int(p_abs)))
            comp = 1 if is_segdup(chrom, int(p_abs)) else 0
            pool = []
            for db in (-1, 0, 1):
                for dr in (-1, 0, 1):
                    k = f'{sb + db}_{rb + dr}_{comp}'
                    if k in z.files:
                        pool.append(z[k])
            if not pool:
                continue
            pool = np.concatenate(pool)
            if len(pool) < 50:
                continue
            p = float((np.searchsorted(pool, age, side='right') + 1) / (len(pool) + 1))
            rho = rho_at(bps, int(p_abs))
            rows.append((int(p_abs), carriers, age, rho, len(pool), p, comp))
        print(f'  {f.name}: {len(rows)} cumulative', flush=True)
    df = pd.DataFrame(rows, columns=['pos', 'carriers', 'age', 'rho', 'null_n', 'p_age', 'segdup'])
    df.to_csv(OUT / f'composite_score_chr{chrom}.csv', index=False)
    print(f'chr{chrom}: scored {len(df)} variants', flush=True)


def region():
    import glob
    from scipy import stats
    frames = []
    for f in sorted(glob.glob(str(OUT / 'composite_score_chr*.csv'))):
        c = int(f.split('chr')[-1].split('.')[0])
        d = pd.read_csv(f); d['chrom'] = c; frames.append(d)
    d = pd.concat(frames, ignore_index=True)
    print(f'scored: {len(d):,} variants | segdup compartment: {100*d.segdup.mean():.1f}%')
    out_rows = []
    for comp in (0, 1):
        sub = d[d.segdup == comp]
        if not len(sub):
            continue
        q05 = sub.p_age.quantile(.05); q01 = sub.p_age.quantile(.01)
        sub = sub.assign(win=(sub.pos // 1_000_000) * 1_000_000)
        g = sub.groupby(['chrom', 'win']).agg(n=('p_age', 'size'),
                                             k5=('p_age', lambda x: (x <= q05).sum()),
                                             k1=('p_age', lambda x: (x <= q01).sum()),
                                             med_rho=('rho', 'median')).reset_index()
        g = g[g.n >= 50].copy()
        g['exp5'] = g.n * .05
        g['enr5'] = g.k5 / g.exp5.clip(lower=.01)
        g['p_val'] = stats.poisson.sf(g.k5 - 1, g.exp5)
        g = g.sort_values('p_val').reset_index(drop=True)
        m = len(g)
        g['q'] = np.minimum.accumulate((g.p_val * m / (g.index + 1))[::-1])[::-1]
        g['segdup'] = comp
        label = 'SEGDUP compartment' if comp else 'NON-SEGDUP compartment'
        print(f'\n=== {label} ({label and len(sub):,} variants, {m} windows, q<0.05: {int((g.q<.05).sum())})')
        print(g.nsmallest(12, 'p_val')[['chrom', 'win', 'n', 'k5', 'exp5', 'enr5', 'p_val', 'q', 'med_rho']].to_string(index=False))
        print('known loci in this compartment:')
        for lab, c, p in (('LCT', 2, 136608646), ('AMY', 1, 104121000), ('SLC45A2', 5, 33951693),
                          ('FADS', 11, 61585000), ('SLC24A5', 15, 48426484)):
            r = g[(g.chrom == c) & (g.win == (p // 1_000_000) * 1_000_000)]
            if len(r):
                r = r.iloc[0]
                rank = int((g.p_val < r.p_val).sum()) + 1
                print(f'  {lab:9s} n={int(r.n)} enr={r.enr5:.2f}x p={r.p_val:.3f} q={r.q:.2f} rank={rank}/{m}')
        out_rows.append(g)
    pd.concat(out_rows, ignore_index=True).to_csv(OUT / 'composite_regional.csv', index=False)
    print('\nwrote composite_regional.csv')


def main():
    mode = sys.argv[1]
    if mode == 'build-null':
        build_null(int(sys.argv[2]) if len(sys.argv) > 2 else 5,
                   int(sys.argv[3]) if len(sys.argv) > 3 else None,
                   sys.argv[4] if len(sys.argv) > 4 else None,
                   int(sys.argv[5]) if len(sys.argv) > 5 else 5)
    elif mode == 'merge-null':
        import glob as _g
        cells = {}
        for f in sorted(_g.glob(str(OUT / 'composite_null_part*.npz'))):
            zz = np.load(f)
            for k in zz.files:
                cells.setdefault(k, []).append(zz[k])
            print('merged', str(f).split('/')[-1], flush=True)
        arr = {k: np.sort(np.concatenate(v)) for k, v in cells.items()}
        np.savez_compressed(OUT / 'composite_null.npz', **arr)
        print(f'wrote composite_null.npz: {len(arr)} cells, {sum(len(v) for v in arr.values()):,} clades')
    elif mode == 'score':
        score(int(sys.argv[2]))
    elif mode == 'region':
        region()


if __name__ == '__main__':
    main()
