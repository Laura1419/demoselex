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
WINDOW = 0          # 0 = whole-chunk null; otherwise +/- bp around the target


def clade_table(ts, chunk_start):
    """(size, age, left, right) for every clade in every tree of this chunk.

    left/right are the tree interval bounds in GENOMIC coordinates (chunk_start added),
    so a null can be scoped to clades from trees overlapping a target position.
    """
    t = ts.tables.nodes.time
    samples = ts.samples()
    sizes, ages, lefts, rights = [], [], [], []
    for tree in ts.trees():
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
        sizes.append(counts[nz]); ages.append(t[nz])
        lefts.append(np.full(len(nz), tree.interval.left + chunk_start))
        rights.append(np.full(len(nz), tree.interval.right + chunk_start))
    return dict(size=np.concatenate(sizes), age=np.concatenate(ages),
                left=np.concatenate(lefts), right=np.concatenate(rights))


def score_all(chrom):
    """Score every variant in every chunk of the chromosome (regional excess input)."""
    rows = []
    chunks = sorted({int(f.name.split('_')[1]) for f in OUT.glob(f'chr{chrom}_*_*.dated.trees')})
    for chunk in chunks:
        f = OUT / f'chr{chrom}_{chunk}_{chunk + 10_000_000}.dated.trees'
        if not f.exists():
            continue
        ts = tskit.load(str(f))
        t = ts.tables.nodes.time
        G = ts.genotype_matrix()
        spos = np.array([s.position for s in ts.sites()]) + chunk
        ct = clade_table(ts, chunk)
        ct['sbin'] = np.round(np.log2(np.clip(ct['size'], 1, None))).astype(int)
        null_by_bin = {b: np.nonzero(ct['sbin'] == b)[0] for b in np.unique(ct['sbin'])}
        mutnode = {}
        for m in ts.mutations():
            p_abs = int(ts.site(m.site).position) + chunk
            if p_abs not in mutnode:
                mutnode[p_abs] = m.node
        for j, p_abs in enumerate(spos):
            carriers = int(G[j].sum())
            if carriers < 3:
                continue
            node = mutnode.get(int(p_abs))
            if node is None:
                continue
            age = float(t[node])
            sbin = int(round(np.log2(carriers)))
            idx = np.concatenate([null_by_bin.get(b, np.array([], dtype=int))
                                  for b in (sbin - 1, sbin, sbin + 1)])
            if len(idx) < 30:
                continue
            pool = ct['age'][idx]
            p = float((np.sum(pool <= age) + 1) / (len(pool) + 1))
            rows.append((int(p_abs), carriers, age, float(np.median(pool)), len(pool), p))
        print(f'  chr{chrom} chunk {chunk}: {len(rows)} cumulative variants', flush=True)
    df = pd.DataFrame(rows, columns=['pos', 'carriers', 'age', 'null_med', 'null_n', 'p_age'])
    out = OUT / f'cladenull_all_chr{chrom}.csv'
    df.to_csv(out, index=False)
    print(f'wrote {out} ({len(df)} variants)')
    return df


def main():
    global WINDOW
    chrom = int(sys.argv[1])
    spec = sys.argv[2]
    if len(sys.argv) > 3:
        WINDOW = int(sys.argv[3])
    if spec == 'all':
        return score_all(chrom)
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
        ct['sbin'] = np.round(np.log2(np.clip(ct['size'], 1, None))).astype(int)
        print(f'  chunk {chunk}: {len(ct["size"])} clades, {len(items)} targets', flush=True)
        null_by_bin = {}
        for b in np.unique(ct['sbin']):
            m = ct['sbin'] == b
            null_by_bin[b] = np.nonzero(m)[0]
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
            idx = np.concatenate([null_by_bin.get(b, np.array([], dtype=int))
                                  for b in (sbin - 1, sbin, sbin + 1)])
            if len(idx) == 0:
                continue
            if WINDOW > 0:          # scope the null to trees overlapping pos +/- WINDOW
                near = (ct['right'][idx] > pos - WINDOW) & (ct['left'][idx] <= pos + WINDOW)
                idx = idx[near]
            if len(idx) < 30:
                continue
            pool = ct['age'][idx]
            p = float((np.sum(pool <= age) + 1) / (len(pool) + 1))
            allrows.append(dict(chrom=chrom, pos=int(pos), carriers=carriers, age=age,
                                null_med=float(np.median(pool)), null_n=len(pool),
                                p_age=p, z_heuristic=zc, freq=fr))
    df = pd.DataFrame(allrows)
    tag = (('top' + spec) if spec.isdigit() else Path(spec).stem) + (f'_w{WINDOW//1000}kb' if WINDOW else '')
    out = OUT / f'cladenull_chr{chrom}_{tag}.csv'
    df.to_csv(out, index=False)
    print(f'wrote {out} ({len(df)} variants)')
    if len(df):
        print('\np_age < 0.05: %d | < 0.01: %d | < 0.001: %d' %
              ((df.p_age < .05).sum(), (df.p_age < .01).sum(), (df.p_age < .001).sum()))
        print(df.nsmallest(12, 'p_age')[['pos', 'carriers', 'age', 'null_med', 'p_age', 'z_heuristic']].to_string(index=False))


if __name__ == '__main__':
    main()
