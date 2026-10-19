#!/usr/bin/env python3
"""omega-style branch-length statistic on dated ARGs (empirical expectation).

Relate's omega compares observed with EXPECTED branch length in recent time bands, the
expectation coming from a coalescent with estimated Ne(t). Here the expectation is taken
EMPIRICALLY from the same ARG: branches are stratified by (time band, log2 descendant-leaf
count), and

    omega_b = sum_bands [ L_{b,j} / E_j(leafbin_b) ] - (#bands with L_{b,j} > 0)

so omega_b > 0 means the branch is longer than comparable branches in recent time - the
sweep signature (a recently-rising allele sits above a clade that coalesces recently, so
its branch is unusually long in the recent bands). Per SNP: omega of the branch carrying
its mutation. This is an empirical analogue of omega (not Relate's exact formula), and is
demography-free because the demography is inside the tree.

Usage: g1k_omega.py <chrom> <array of positions> [<lo> <hi>]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
POPFILE = OUT / 'samples_population.txt'
BANDS = np.array([0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                  10_000, 20_000, 50_000, 100_000, 1_000_000], dtype=float)


def eur_haps():
    sel = []
    with open(POPFILE) as f:
        next(f)
        for i, line in enumerate(f):
            p = line.split()
            if len(p) >= 7 and p[6] == 'EUR':
                sel.append(i)
    return [h for i in sel for h in (2 * i, 2 * i + 1)]


def chunk_for(chrom, pos):
    for f in sorted(OUT.glob(f'chr{chrom}_*.dated.trees')):
        s = int(f.name.split('_')[1]); e = int(f.name.split('_')[2].split('.')[0])
        if s <= pos < e:
            return f, s
    raise FileNotFoundError((chrom, pos))


def leaf_counts(tree, t, n_nodes, samples):
    cnt = np.zeros(n_nodes, dtype=np.int32)
    cnt[samples] = 1
    par = tree.parent_array
    for v in np.argsort(t, kind='stable'):
        p = par[v]
        if p != -1:
            cnt[p] += cnt[v]
    return cnt


def compute(chrom, lo, hi):
    haps = eur_haps()
    f, start = chunk_for(chrom, lo)
    ts = tskit.load(str(f)).simplify(samples=haps)
    t = ts.tables.nodes.time
    samples = ts.samples()
    rows = []     # (branch_id, band, leafbin, length)
    muts = []     # (pos, branch_id)
    bid = 0
    branch_leafbin = []
    for ti, tree in enumerate(ts.trees()):
        left, right = tree.interval.left + start, tree.interval.right + start
        if right <= lo or left >= hi:
            continue
        cnt = leaf_counts(tree, t, ts.num_nodes, samples)
        par = tree.parent_array
        edge_of_child = {}
        for c in np.where(par != -1)[0]:
            p = par[c]
            tc, tp = t[c], t[p]
            if tp <= 0 or not np.isfinite(tp):
                continue
            lb = int(np.floor(np.log2(max(int(cnt[c]), 1))))
            branch_leafbin.append(lb)
            edge_of_child[int(c)] = bid
            for bi in range(len(BANDS) - 1):
                a, b = BANDS[bi], BANDS[bi + 1]
                ov = min(tp, b) - max(tc, a)
                if ov > 0:
                    rows.append((bid, bi, lb, float(ov)))
            bid += 1
        for m in tree.mutations():
            pp = ts.site(m.site).position + start
            if lo <= pp < hi:
                b = edge_of_child.get(int(m.node))
                if b is not None:
                    muts.append((int(pp), b))
    df = pd.DataFrame(rows, columns=['bid', 'band', 'leafbin', 'len'])
    E = df.groupby(['band', 'leafbin']).len.agg(['mean', 'size'])
    E = E[E['size'] >= 20]['mean'].to_dict()
    df['E'] = [E.get((b, l), np.nan) for b, l in zip(df.band, df.leafbin)]
    df = df.dropna(subset=['E'])
    df['ratio'] = df.len / df.E
    per_branch = df.groupby('bid').agg(omega=('ratio', lambda x: x.sum() - len(x)))
    om = per_branch.omega.to_dict()
    n_bands = df.groupby('bid').size().to_dict()
    out = pd.DataFrame([dict(pos=p, omega=om.get(b, np.nan), bands=n_bands.get(b, 0)) for p, b in muts])
    return out, len(df), len(E)


def main():
    chrom = int(sys.argv[1])
    targets = [(x.split(':')[0], int(x.split(':')[1])) for x in sys.argv[2].split(',')]
    print(f'chr{chrom} | empirical omega (branch length vs comparable branches, recent bands)')
    rows = []
    # group targets by 10 Mb chunk to reuse the loaded ARG
    groups = {}
    for lab, p in targets:
        groups.setdefault(int(p // 10_000_000), []).append((lab, p))
    for chunk, items in sorted(groups.items()):
        lo = chunk; hi = chunk + 10_000_000
        out, nb, ns = compute(chrom, lo, hi)
        if not len(out):
            print(f'  chunk {chunk}: no mutations'); continue
        for lab, p in items:
            m = out[np.abs(out.pos - p) <= 1]
            if not len(m):
                m = out[np.abs(out.pos - p) <= 1000]
            if not len(m):
                rows.append(dict(locus=lab, pos=p, status='absent')); continue
            r = m.iloc[np.abs(m.pos - p).argmin()]
            pct = 100 * (out.omega < r.omega).mean()
            rows.append(dict(locus=lab, pos=int(r.pos), omega=round(float(r.omega), 3),
                             percentile=round(float(pct), 1), window_median=round(float(out.omega.median()), 3),
                             n_mut=len(out), strata=ns))
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == '__main__':
    main()
