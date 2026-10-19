#!/usr/bin/env python3
"""Local, ARG-derived backgrounds (no global population labels).

For a locus, each haplotype's **local background** is defined from the local tree as the
most recent clade containing it whose root is younger than tau:

    root(s) = deepest ancestor of sample s with time <= tau  (else s itself)

So the partition is learned per locus from the genealogy, at a chosen time depth, rather
than imposed from external population labels.

Modes
  locus <chrom> <pos> [tau]     partition + per-cluster allele frequency at one site
  window <chrom> <lo> <hi> [tau] partition at several positions, cluster sizes, and
                                cross-tabulation against the global labels
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
POPFILE = OUT / 'samples_population.txt'
TAU = 200.0


def labels():
    pops, sups = [], []
    with open(POPFILE) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) >= 7:
                pops.append(p[5]); sups.append(p[6])
    return pops, sups


def chunk_for(chrom, pos):
    start = int(pos // 10_000_000) * 10_000_000
    for f in OUT.glob(f'chr{chrom}_{start}_*.dated.trees'):
        return f, start
    # fall back to any chunk containing the position
    for f in sorted(OUT.glob(f'chr{chrom}_*.dated.trees')):
        s = int(f.name.split('_')[1]); e = int(f.name.split('_')[2].split('.')[0])
        if s <= pos < e:
            return f, s
    raise FileNotFoundError(f'no chunk containing chr{chrom}:{pos}')


def local_roots(ts, pos, tau):
    """Cluster root node per haplotype at `pos` (most recent ancestor younger than tau)."""
    t = ts.tables.nodes.time
    tree = ts.at(pos)
    par = tree.parent_array
    out = np.zeros(ts.num_samples, dtype=np.int64)
    for k, s in enumerate(ts.samples()):
        v = s
        while par[v] != -1 and t[par[v]] <= tau:
            v = par[v]
        out[k] = v
    return out


def summarize(roots, ts, pos, tau, allele=None):
    uniq, inv, counts = np.unique(roots, return_inverse=True, return_counts=True)
    big = [(int(u), int(c)) for u, c in zip(uniq, counts) if c >= 20]
    pops, sups = labels()
    # cross-tab: for the largest local clusters, what global superpopulation are they?
    rows = []
    for u, c in sorted(big, key=lambda x: -x[1])[:10]:
        idx = np.where(inv == np.where(uniq == u)[0][0])[0]
        samp = idx // 2
        sup = pd.Series([sups[i] for i in samp]).value_counts().head(3).to_dict()
        pop = pd.Series([pops[i] for i in samp]).value_counts().head(3).to_dict()
        rows.append(dict(cluster=int(u), n_hap=int(c), n_samples=int(len(set(samp))),
                         superpops=sup, pops=pop))
    return rows, len(big), len(uniq)


def main():
    mode = sys.argv[1]
    chrom = int(sys.argv[2])
    if mode == 'locus':
        pos = int(sys.argv[3])
        tau = float(sys.argv[4]) if len(sys.argv) > 4 else TAU
        f, start = chunk_for(chrom, pos)
        ts = tskit.load(str(f))
        rel = pos - start
        roots = local_roots(ts, rel, tau)
        j = int(np.argmin(np.abs(np.array([s.position for s in ts.sites()]) - rel)))
        G = ts.genotype_matrix()[j]
        print(f'chr{chrom}:{pos} tau={tau} gen | site {int(ts.site(j).position)+start} '
              f'| carriers {int(G.sum())}/{len(G)} ({G.mean():.3f})')
        rows, nbig, ntot = summarize(roots, ts, rel, tau)
        print(f'local clusters (>=20 haps): {nbig} of {ntot}')
        df = pd.DataFrame(rows)
        # per-cluster allele frequency
        freq = []
        for r in rows:
            m = roots == r['cluster']
            freq.append(round(float(G[m].mean()), 3))
        df['allele_freq'] = freq
        print(df.to_string(index=False))
        print('\ncarrier distribution across local clusters:',
              {int(u): int(G[roots == u].sum()) for u, _ in
               sorted([(int(u), int(c)) for u, c in zip(*np.unique(roots, return_counts=True))
                       if c >= 20], key=lambda x: -x[1])[:8]})
    elif mode == 'confine':
        tau = float(sys.argv[3]) if len(sys.argv) > 3 else TAU
        loci = sys.argv[4].split(',')
        print(f'tau={tau} gen | background confinement of carriers (label-free)')
        rows = []
        for spec in loci:
            lab, p_ = spec.split(':')
            pos = int(p_)
            f, start = chunk_for(chrom, pos)
            ts = tskit.load(str(f))
            rel = pos - start
            roots = local_roots(ts, rel, tau)
            j = int(np.argmin(np.abs(np.array([s.position for s in ts.sites()]) - rel)))
            G = ts.genotype_matrix()[j]
            u, cnt = np.unique(roots, return_counts=True)
            car = {int(x): int(G[roots == x].sum()) for x in u}
            tot = int(G.sum())
            # carriers' cluster distribution
            top = sorted(car.items(), key=lambda kv: -kv[1])[:5]
            n_clusters_carrying = sum(1 for k, v in car.items() if v >= max(1, 0.05 * tot))
            big = {int(x): int(c) for x, c in zip(u, cnt) if c >= 20}
            dom = top[0][0] if top else None
            dom_freq = (car[dom] / big[dom]) if dom in big and big[dom] else float('nan')
            rows.append(dict(locus=lab, pos=int(ts.site(j).position) + start, carriers=tot,
                             freq_global=round(tot / len(G), 3),
                             clusters_20plus=len(big),
                             clusters_ge5pct_carriers=n_clusters_carrying,
                             dominant_share=round(top[0][1] / tot, 3) if tot else float('nan'),
                             dom_cluster_freq=round(dom_freq, 3),
                             dom_cluster_n=big.get(dom, 0)))
        print(pd.DataFrame(rows).to_string(index=False))
    elif mode == 'window':
        lo, hi = int(sys.argv[3]), int(sys.argv[4])
        tau = float(sys.argv[5]) if len(sys.argv) > 5 else TAU
        f, start = chunk_for(chrom, lo)
        ts = tskit.load(str(f))
        print(f'window chr{chrom}:{lo}-{hi} tau={tau} gen')
        for p in range(lo, hi, 100_000):
            roots = local_roots(ts, p - start, tau)
            rows, nbig, ntot = summarize(roots, ts, p - start, tau)
            print(f'  {p}: {nbig} clusters >=20 haps (of {ntot}); largest '
                  f'{rows[0]["n_hap"] if rows else 0} haps, superpops {rows[0]["superpops"] if rows else {}}')


if __name__ == '__main__':
    main()
