#!/usr/bin/env python3
"""Population-restricted clade-age statistics from the 1000G dated ARGs.

The global analysis dilutes population-specific signal: carrier clades span
populations (deeper MRCA), frequencies are diluted, and the frequency-matched null
fills with young population-specific variants. This restricts the genealogy to a
sample subset (by super-population) via tskit simplify, then recomputes clade ages.

Usage: g1k_pop.py stats <chrom> <SUPERPOP|comma-list-of-pops>
       g1k_pop.py check <chrom> <SUPERPOP> <label:pos[,label:pos...]>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
POPFILE = OUT / 'samples_population.txt'
HALF = 250_000


def pop_index(pop_spec):
    """Sample indices (0-based, matching VCF order) for a superpopulation or pop list."""
    table = []
    with open(POPFILE) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 7:
                table.append((p[1], p[5], p[6]))       # sample, pop, superpop
    wanted = set(pop_spec.split(','))
    return [i for i, (s, pop, sup) in enumerate(table) if pop in wanted or sup in wanted]


def stats(chrom, pop_spec):
    import tskit, glob
    sel = pop_index(pop_spec)
    haps = [h for i in sel for h in (2 * i, 2 * i + 1)]
    tag = pop_spec.replace(',', '_')
    rows = []
    for f in sorted(glob.glob(str(OUT / f'chr{chrom}_*.dated.trees'))):
        start = int(f.split('_')[-2])
        ts = tskit.load(f)
        sub = ts.simplify(samples=haps)
        t = sub.tables.nodes.time
        G = sub.genotype_matrix()
        spos = np.array([s.position for s in sub.sites()]) + start
        bps = np.array(list(sub.breakpoints())) + start
        node = {}
        for m in sub.mutations():
            k = int(sub.site(m.site).position) + start
            node.setdefault(k, m.node)
        for j, p in enumerate(spos):
            n = int(G[j].sum())
            if n < 3 or int(p) not in node:
                continue
            age = float(t[node[int(p)]])
            lo = np.searchsorted(bps, p - HALF); hi = np.searchsorted(bps, p + HALF)
            rho = (hi - lo) / (2 * HALF / 1000.0)
            rows.append((int(p), n, age, rho))
        print(f'  {Path(f).name}: {len(sub.num_sites * [0])} sites, {len(rows)} rows', flush=True)
    df = pd.DataFrame(rows, columns=['pos', 'carriers', 'age', 'rho'])
    out = OUT / f'chr{chrom}_{tag}_stats.csv'
    df.to_csv(out, index=False)
    print(f'wrote {out} ({len(df)} variants; {len(sel)} samples, {len(haps)} haplotypes)', flush=True)
    return df


def check(chrom, pop_spec, targets):
    import tskit, glob
    tag = pop_spec.replace(',', '_')
    f = OUT / f'chr{chrom}_{tag}_stats.csv'
    d = pd.read_csv(f) if f.exists() else stats(chrom, pop_spec)
    n_h = int(d.carriers.max())  # not the true n; recompute from samples
    sel = pop_index(pop_spec)
    n_h = 2 * len(sel)
    d['freq'] = d.carriers / n_h
    d['fb'] = pd.qcut(d.freq, 12, duplicates='drop')
    cur = d.groupby('fb', observed=True).age.agg(['median', 'std', 'count'])
    q5 = d.age.quantile(.05)
    print(f'chr{chrom} {pop_spec}: {len(d):,} variants | {n_h} haplotypes | median age {d.age.median():.0f} gen')
    for spec in targets:
        lab, pos = spec.split(':')
        pos = int(pos)
        r = d[np.abs(d.pos - pos) <= 1]
        if not len(r):
            r = d[np.abs(d.pos - pos) <= 1000]
        if not len(r):
            print(f'  {lab}: absent'); continue
        r = r.iloc[np.abs(r.pos - pos).argmin()]
        b = pd.cut([r.freq], d.fb.cat.categories)[0]
        sub = d[d.fb == b]
        pct = 100 * (sub.age < r.age).mean()
        w = d[(d.pos >= pos - 1_000_000) & (d.pos <= pos + 1_000_000)]
        k = int((w.age <= q5).sum()); exp = len(w) * 0.05
        from scipy import stats as st
        pv = st.poisson.sf(k - 1, exp) if exp > 0 else np.nan
        print(f'  {lab:22s} pos={int(r.pos)} freq={r.freq:.3f} n={int(r.carriers):5d} age={r.age:7.0f} gen '
              f'| pct={pct:5.1f}% | window {k}/{exp:.0f} = {k/exp if exp else np.nan:.1f}x p={pv:.1e}')


def main():
    mode = sys.argv[1]
    chrom = int(sys.argv[2])
    pop = sys.argv[3]
    if mode == 'stats':
        stats(chrom, pop)
    elif mode == 'check':
        check(chrom, pop, sys.argv[4].split(','))


if __name__ == '__main__':
    main()
