#!/usr/bin/env python3
"""Clade-size-matched age null (corrected within-tree null).

A mutation in a tsinfer ARG is placed ON a clade, so "carriers form a clade" is a
tautology and cannot be tested. What IS testable: is this clade unusually YOUNG for
how many samples it contains? The local genealogy supplies its own null: enumerate
every clade in the tree (every node's descendant sample set), and compare the observed
clade's age to the distribution of ages of clades with the same size.

This is demography-free (the demography is in the tree), needs no simulation, and is
the natural ARG analogue of an excess-of-recent-coalescence statistic.

Usage: arg_cladenull.py <chrom> <N|targets.csv> [window_kb]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/arggenome'


def clade_table(ts, chunk_start):
    """(size, age, midpoint) for every clade present in every tree of this chunk."""
    t = ts.tables.nodes.time
    samples = ts.samples()
    sset = np.zeros(ts.num_nodes, dtype=bool)
    sset[samples] = True
    rows = []
    for tree in ts.trees():
        par = tree.parent_array
        counts = np.zeros(ts.num_nodes, dtype=np.int32)
        for s in samples:
            v = s
            while v != -1:
                counts[v] += 1
                v = par[v]
        mid = (tree.interval.left + tree.interval.right) / 2 + chunk_start
        for v in np.nonzero(counts)[0]:
            rows.append((int(counts[v]), float(t[v]), float(mid)))
    return pd.DataFrame(rows, columns=['size', 'age', 'mid'])


def main():
    chrom = int(sys.argv[1])
    spec = sys.argv[2]
    z = pd.read_csv(OUT / 'genomewide_z.csv')
    if spec.isdigit():
        sub = z[z.chrom == chrom].nsmallest(int(spec), 'z')
    else:
        sub = pd.read_csv(spec)
    sub = sub[sub.chrom == chrom] if 'chrom' in sub.columns else sub
    by_chunk = {}
    for r in sub.itertuples():
        by_chunk.setdefault(int(r.chunk), []).append((int(r.pos), float(r.z), float(r.freq)))
    allrows = []
    for chunk, items in sorted(by_chunk.items()):
        f = OUT / f'chr{chrom}_{chunk}_{chunk + 10_000_000}.dated.trees'
        if not f.exists():
            continue
        ts = tskit.load(str(f))
        t = ts.tables.nodes.time
        G = ts.genotype_matrix()
        spos = np.array([s.position for s in ts.sites()])
        ct = clade_table(ts, chunk)
        print(f'  chunk {chunk}: {len(ct)} clades, {len(items)} targets', flush=True)
        # size-stratified null (log2 size bins, smoothed over neighbours)
        ct['sbin'] = np.round(np.log2(ct['size'].clip(1))).astype(int)
        null_by_bin = {b: g.age.values for b, g in ct.groupby('sbin')}
        for pos, zc, fr in items:
            rel = pos - chunk
            j = int(np.argmin(np.abs(spos - rel)))
            if abs(spos[j] - rel) > 1:
                continue
            carriers = int(G[j].sum())
            if carriers < 3:
                continue
            mut = [m for m in ts.mutations() if m.site == j]
            if not mut:
                continue
            node = mut[0].node
            age = float(t[node])
            sbin = int(round(np.log2(max(carriers, 1))))
            pool = np.concatenate([null_by_bin.get(b, np.array([])) for b in (sbin - 1, sbin, sbin + 1)])
            if len(pool) < 30:
                continue
            p = float((np.sum(pool <= age) + 1) / (len(pool) + 1))
            allrows.append(dict(chrom=chrom, pos=int(pos), carriers=carriers, age=age,
                                null_med=float(np.median(pool)), null_n=len(pool),
                                p_age=p, z_heuristic=zc, freq=fr))
    df = pd.DataFrame(allrows)
    tag = ('top' + spec) if spec.isdigit() else Path(spec).stem
    out = OUT / f'cladenull_chr{chrom}_{tag}.csv'
    df.to_csv(out, index=False)
    print(f'wrote {out} ({len(df)} variants)')
    if len(df):
        print('\np_age < 0.05: %d | < 0.01: %d | < 0.001: %d' %
              ((df.p_age < .05).sum(), (df.p_age < .01).sum(), (df.p_age < .001).sum()))
        print(df.nsmallest(12, 'p_age')[['pos', 'carriers', 'age', 'null_med', 'p_age', 'z_heuristic']].to_string(index=False))


if __name__ == '__main__':
    main()
