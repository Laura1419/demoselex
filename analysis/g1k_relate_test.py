#!/usr/bin/env python3
"""Relate-style selection test, implemented empirically.

CORRECTION OF RECORD: Relate (Speidel et al. 2019) has NO omega statistic. Its selection
test conditions on k0 = the number of lineages present when the mutation arises, and uses
n0 = the number of present-day carriers as the statistic, with an analytic null robust to
population-size change (p_R). omega/OmegaPlus is a different, haplotype-based method.

Implementation here: for every mutation,
    n0 = number of descendant samples of the mutation node  (carriers)
    k0 = 1 + number of internal nodes in the local tree with time > t_node  (lineages at origin)
Null (empirical, genome/chromosome-wide, replacing Relate's analytic null):
    the distribution of n0 among mutations with the same k0.
    p_R = fraction of that null with n0' >= n0   (one-sided: spread faster than expected)

Usage: g1k_relate_test.py <chrom> <lo> <hi> [label:pos,...]
       g1k_relate_test.py null <chrom> [stride]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
POPFILE = OUT / 'samples_population.txt'


def eur_haps():
    sel = []
    with open(POPFILE) as f:
        next(f)
        for i, line in enumerate(f):
            p = line.split()
            if len(p) >= 7 and p[6] == 'EUR':
                sel.append(i)
    return [h for h in (h for i in sel for h in (2 * i, 2 * i + 1))]


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


def collect(chrom, lo=None, hi=None, stride=1):
    """(pos, n0, k0) for every mutation in the range."""
    rows = []
    for f in chunks_for(chrom):
        start = int(f.name.split('_')[1]); end = int(f.name.split('_')[2].split('.')[0])
        if lo is not None and (end <= lo or start >= hi):
            continue
        ts = tskit.load(str(f)).simplify(samples=eur_haps())
        t = ts.tables.nodes.time
        samples = ts.samples()
        for ti, tree in enumerate(ts.trees()):
            if stride > 1 and ti % stride:
                continue
            l_, r_ = tree.interval.left + start, tree.interval.right + start
            if lo is not None and (r_ <= lo or l_ >= hi):
                continue
            cnt = leaf_counts(tree, t, ts.num_nodes, samples)
            par = tree.parent_array
            internal_times = np.sort(t[[v for v in np.where(par != -1)[0] if (par[v] != -1 or True)
                                        and t[v] > 0 and t[v] < np.inf]]) if False else None
            # internal nodes of this tree: nodes with >=1 child and a finite time
            kids = np.where(par != -1)[0]
            itimes = np.sort(np.array([t[v] for v in np.unique(par[kids])]))
            nn = len(itimes)
            for m in tree.mutations():
                pp = ts.site(m.site).position + start
                if lo is not None and not (lo <= pp < hi):
                    continue
                tm = t[m.node]
                k0 = 1 + int(nn - np.searchsorted(itimes, tm, side='right'))
                rows.append((int(pp), int(cnt[m.node]), k0))
    df = pd.DataFrame(rows, columns=['pos', 'n0', 'k0'])
    # a site can carry several mutations (recurrences); keep the one with most carriers
    return df.sort_values('n0').drop_duplicates('pos', keep='last').sort_values('pos')


def build_null(chrom, stride=20):
    df = collect(chrom, stride=stride)
    df.to_csv(OUT / f'relate_null_chr{chrom}.csv', index=False)
    print(f'chr{chrom}: null from {len(df):,} mutations; k0 range {df.k0.min()}-{df.k0.max()}')


def scan(chrom, lo, hi, targets):
    nulls = []
    for f in sorted(OUT.glob('relate_null_chr*.csv')):
        nulls.append(pd.read_csv(f))
    nul = pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame()
    # k0 strata: exact k0, with neighbouring fallback
    byk = {int(k): g.n0.values for k, g in nul.groupby('k0')} if len(nul) else {}
    df = collect(chrom, lo, hi)
    if not len(df):
        print('no mutations'); return
    rows = []
    for r in df.itertuples():
        pool = np.concatenate([byk.get(k, np.array([])) for k in (r.k0 - 1, r.k0, r.k0 + 1)])
        if len(pool) < 50:
            rows.append(dict(pos=r.pos, n0=r.n0, k0=r.k0, pR=np.nan, null_n=len(pool)))
            continue
        pR = float((np.sum(pool >= r.n0) + 1) / (len(pool) + 1))
        rows.append(dict(pos=r.pos, n0=r.n0, k0=r.k0, pR=pR, null_n=len(pool)))
    out = pd.DataFrame(rows)
    if not len(out):
        print('no testable mutations'); return
    # analytic null: given k0 lineages, the carriers' sample fraction is Beta(1, k0-1) so
    # P(fraction >= x) = (1-x)^(k0-1). No resolution limit.
    n_sample = 1266
    out['x'] = out.n0 / n_sample
    out['pR_analytic'] = np.power(np.clip(1 - out.x, 1e-12, 1), np.maximum(out.k0 - 1, 1))
    out['neglog10p_analytic'] = -np.log10(out.pR_analytic.clip(lower=1e-300))
    out['neglog10p'] = -np.log10(out.pR.clip(lower=1e-300))
    print(f'chr{chrom}:{lo}-{hi} | {len(out):,} tested | median pR_emp {out.pR.median():.3f}'
          f' | frac pR_emp<0.01 {100*(out.pR<0.01).mean():.2f}%'
          f' | frac pR_analytic<1e-4 {100*(out.pR_analytic<1e-4).mean():.2f}%')
    for lab, p in targets:
        m = out[np.abs(out.pos - p) <= 1]
        if not len(m):
            m = out[np.abs(out.pos - p) <= 1000]
        if not len(m):
            print(f'  {lab}: absent'); continue
        r = m.iloc[np.abs(m.pos - p).argmin()]
        print(f'  {lab:10s} pos={int(r.pos)} n0={int(r.n0)} k0={int(r.k0)} '
              f'pR_emp={r.pR:.2e} | pR_analytic={r.pR_analytic:.2e} '
              f'(-log10p {r.neglog10p_analytic:.1f})')
    return out


def main():
    mode = sys.argv[1]
    if mode == 'null':
        build_null(int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20)
    else:
        chrom = int(sys.argv[1]); lo = int(sys.argv[2]); hi = int(sys.argv[3])
        targets = [(x.split(':')[0], int(x.split(':')[1])) for x in sys.argv[4].split(',')] if len(sys.argv) > 4 else []
        scan(chrom, lo, hi, targets)


if __name__ == '__main__':
    main()
