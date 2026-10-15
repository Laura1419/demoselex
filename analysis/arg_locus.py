#!/usr/bin/env python3
"""ARG-based allele-age test vs frequency (per chromosome).

For each site in the dated ARG, the derived allele's age is the tsdate time of the
node the mutation sits on (the MRCA of its carriers for a unique mutation). Under
neutrality allele age rises with frequency; a recent sweep produces a derivative
allele that is COMMON but YOUNG, i.e. below the neutral age-vs-frequency curve.

Null: the empirical age-vs-frequency relation on the same chromosome, in frequency
bins. Output: per-site table, plus the rank/percentile of target-region sites within
their frequency bin.

Usage: arg_locus.py <dataset> <chrom> <label:start-end[,label:start-end...]>
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

OUT = Path(__file__).resolve().parents[1] / 'data/derived/arg'
NBINS = 12
MIN_SAMPLES_PER_BIN = 50


def site_table(dataset, chrom):
    ts = tskit.load(str(OUT / f'{dataset}_chr{chrom}.dated.trees'))
    time = ts.tables.nodes.time
    mut = ts.tables.mutations
    pos = np.array([s.position for s in ts.sites()])
    # derived allele counts per site
    G = ts.genotype_matrix()                     # (sites, samples) 0/1
    n = G.shape[1]
    cnt = G.sum(1)
    freq = cnt / n
    # mutation age per site: min over mutations at that site (recurrent -> youngest)
    age = np.full(ts.num_sites, np.nan)
    for s, nd in zip(mut.site, mut.node):
        t = time[nd]
        if np.isnan(age[s]) or t < age[s]:
            age[s] = t
    return pd.DataFrame(dict(pos=pos, freq=freq, age=age, n_carrier=cnt))


def null_curve(df):
    d = df[(df.freq > 0) & (df.age > 0)].copy()
    d['bin'] = pd.qcut(d.freq, NBINS, duplicates='drop')
    g = d.groupby('bin', observed=True).age.agg(['median', 'count', 'std'])
    return d, g


def main():
    dataset, chrom = sys.argv[1], int(sys.argv[2])
    targets = []
    for spec in sys.argv[3].split(','):
        lab, rng = spec.split(':')
        lo, hi = map(int, rng.split('-'))
        targets.append((lab, lo, hi))
    df = site_table(dataset, chrom)
    d, curve = null_curve(df)
    print('chromosome-wide age-vs-frequency bins (null):')
    print(curve.to_string())
    rows = []
    for lab, lo, hi in targets:
        sub = d[(d.pos >= lo) & (d.pos <= hi)]
        for r in sub.itertuples():
            b = r.bin
            med = curve.loc[b, 'median'] if b in curve.index else np.nan
            sd = curve.loc[b, 'std'] if b in curve.index else np.nan
            z = (r.age - med) / sd if sd and sd > 0 else np.nan
            rows.append(dict(label=lab, pos=int(r.pos), freq=round(float(r.freq), 4),
                             age_gen=round(float(r.age), 1), n_carrier=int(r.n_carrier),
                             bin_median=round(float(med), 1) if np.isfinite(med) else np.nan,
                             z_vs_bin=round(float(z), 2) if np.isfinite(z) else np.nan))
    out = pd.DataFrame(rows)
    print('\ntarget-region sites:')
    print(out.to_string(index=False) if len(out) else '(no sites in target regions)')
    out.to_csv(OUT / f'{dataset}_chr{chrom}_locus.csv', index=False)
    summary = dict(dataset=dataset, chrom=chrom, n_sites=int(len(df)),
                   n_target=int(len(out)),
                   n_young_z_lt_m2=int((out.z_vs_bin < -2).sum()) if len(out) else 0,
                   min_z=float(out.z_vs_bin.min()) if len(out) else None)
    json.dump(summary, open(OUT / f'{dataset}_chr{chrom}_locus_summary.json', 'w'), indent=1)
    print(json.dumps(summary))


if __name__ == '__main__':
    main()
