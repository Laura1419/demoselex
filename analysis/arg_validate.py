#!/usr/bin/env python3
"""Run the modern-ARG clade-age test across a list of windows and report, for each
target variant, its frequency, carrier-clade age, and its rank among common variants.

Usage: arg_validate.py <chrom> <lo> <hi> <label:pos> [<label:pos> ...]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
OUT = ROOT / 'data/derived/argmodern'
COMMON_FREQ = 0.30


def report(chrom, lo, hi, targets):
    ts = tskit.load(str(OUT / f'chr{chrom}_{lo}_{hi}.dated.trees'))
    t = ts.tables.nodes.time
    mut = ts.tables.mutations
    pos = np.array([s.position for s in ts.sites()]) + lo
    G = ts.genotype_matrix()
    freq = G.sum(1) / G.shape[1]
    age = np.full(ts.num_sites, np.nan)
    for s, nd in zip(mut.site, mut.node):
        v = t[nd]
        if not np.isfinite(age[s]) or v < age[s]:
            age[s] = v
    df = pd.DataFrame(dict(pos=pos, freq=freq, age=age))
    d = df[(df.freq > 0.02) & (df.age > 0)]
    comp = d[d.freq >= COMMON_FREQ]
    rows = []
    for lab, p in targets:
        m = df[np.abs(df.pos - p) <= 1]
        if not len(m):
            rows.append(dict(target=lab, pos=p, status='ABSENT'))
            continue
        r = m.iloc[0]
        rank = int((comp.age < r.age).sum()) + 1
        rows.append(dict(target=lab, pos=int(r.pos), freq=round(float(r.freq), 3),
                         age_gen=round(float(r.age), 1), age_yr=int(r.age * 29),
                         rank_among_common=rank, n_common=len(comp),
                         common_median_age=round(float(comp.age.median()), 1)))
    out = pd.DataFrame(rows)
    print(f'--- chr{chrom}:{lo}-{hi}  (ARG sites {ts.num_sites}, samples {ts.num_samples}, common n={len(comp)})')
    print(out.to_string(index=False) if len(out) else '(no sites in window)')
    return rows


def main():
    chrom, lo, hi = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    targets = [(s.split(':')[0], int(s.split(':')[1])) for s in sys.argv[4:]]
    report(chrom, lo, hi, targets)


if __name__ == '__main__':
    main()
