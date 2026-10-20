#!/usr/bin/env python3
"""Scan a BED of loci (gene names in column 4) against the ARG-derived metrics.

For each locus: how its variants sit in the young-clade distribution (population-restricted),
and optionally the local-background breadth of its most extreme variant.

Inputs
  BED: chrom, start, end, name[, ...]      (0-based half-open, GRCh38, as UCSC/cosigt use)
  scores: data/derived/g1k/composite_score_chr<N>.csv  (per-variant p_age, EUR)
  ARG:    data/derived/g1k/chr<N>_<s>_<e>.dated.trees  (for breadth, --breadth)

Usage
  g1k_bed_scan.py loci.bed [--pop EUR] [--breadth] [--tau 200] [--out loci_scan.csv]
"""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
REPO = ROOT
SEGDUP = ROOT / 'data/raw/annotation/segdup/genomicSuperDups_hg38.txt.gz'
POPFILE = OUT / 'samples_population.txt'

_MASK = None


def segdup_mask():
    global _MASK
    if _MASK is not None:
        return _MASK
    from collections import defaultdict
    iv = defaultdict(list)
    with gzip.open(SEGDUP, 'rt') as f:
        for line in f:
            p = line.split('\t')
            if len(p) < 4 or not p[1].startswith('chr'):
                continue
            c = p[1].replace('chr', '')
            if c.isdigit() and 1 <= int(c) <= 22:
                iv[int(c)].append((int(p[2]), int(p[3])))
    m = {}
    for c, lst in iv.items():
        lst.sort(); merged = []
        for a, b in lst:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        m[c] = (np.array([x[0] for x in merged]), np.array([x[1] for x in merged]))
    _MASK = m
    return m


def is_segdup(chrom, pos):
    m = segdup_mask().get(chrom)
    if m is None:
        return False
    s, e = m
    i = np.searchsorted(s, pos, side='right') - 1
    return bool(i >= 0 and pos <= e[i])


def read_bed(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            p = line.split('\t')
            if len(p) < 3:
                continue
            c = p[0].replace('chr', '')
            if not c.isdigit():
                continue
            rows.append(dict(chrom=int(c), start=int(p[1]), end=int(p[2]),
                             name=(p[3].strip() if len(p) > 3 else f'locus_{len(rows)}')))
    return pd.DataFrame(rows)


def load_scores(chrom, pop='EUR'):
    f = OUT / f'composite_score_chr{chrom}.csv'
    if not f.exists():
        return None
    return pd.read_csv(f)


def score_region(df, scores):
    """Region-level summary from precomputed per-variant p_age (bottom-5% excess)."""
    if scores is None or not len(scores) or 'p_age' not in scores:
        return dict(n_variants=0)
    sub = scores[(scores.pos >= df.start) & (scores.pos < df.end)]
    if not len(sub):
        return dict(n_variants=0)
    q5 = scores.p_age.quantile(.05)
    k = int((sub.p_age <= q5).sum())
    exp = len(sub) * 0.05
    from scipy import stats as st
    pv = st.poisson.sf(k - 1, exp) if exp > 0 else np.nan
    best = sub.nsmallest(1, 'p_age').iloc[0]
    return dict(n_variants=len(sub), k5=k, expected=round(exp, 1),
                enrichment=round(k / exp, 2) if exp else np.nan,
                poisson_p=float(pv) if pv == pv else np.nan,
                best_pos=int(best.pos), best_p_age=round(float(best.p_age), 5),
                best_age=round(float(best.age), 0), best_pct=round(100 * float(best.p_age), 2))


def pop_haps(pop='EUR'):
    sel = []
    with open(POPFILE) as f:
        next(f)
        for i, line in enumerate(f):
            p = line.split()
            if len(p) >= 7 and p[6] == pop:
                sel.append(i)
    return [h for i in sel for h in (2 * i, 2 * i + 1)]


def breadth_at(chrom, pos, tau=200.0, pop='EUR'):
    """Local-background breadth of the carriers of the mutation nearest pos (EUR)."""
    import tskit
    f = None
    for g in sorted(OUT.glob(f'chr{chrom}_*.dated.trees')):
        s = int(g.name.split('_')[1]); e = int(g.name.split('_')[2].split('.')[0])
        if s <= pos < e:
            f, start = g, s
            break
    if f is None:
        return dict(breadth=np.nan)
    ts = tskit.load(str(f))
    if pop:
        ts = ts.simplify(samples=pop_haps(pop))
    t = ts.tables.nodes.time
    rel = pos - start
    j = int(np.argmin(np.abs(np.array([s.position for s in ts.sites()]) - rel)))
    site = ts.site(j)
    muts = [m for m in ts.mutations() if m.site == j]
    if not muts:
        return dict(breadth=np.nan)
    tree = ts.at(site.position)
    par = tree.parent_array
    samples = ts.samples()
    # carriers = samples below m.node
    children = {}
    for v in np.where(par != -1)[0]:
        children.setdefault(int(par[v]), []).append(int(v))

    def sub_size(u):
        tot = 0
        st = [int(u)]
        while st:
            x = st.pop(); tot += 1
            st.extend(children.get(x, []))
        return tot
    m = max(muts, key=lambda x: sub_size(x.node))
    below = set()
    stack = [int(m.node)]
    while stack:
        u = stack.pop(); below.add(u)
        stack.extend(children.get(u, []))
    carriers = [s for s in samples if s in below]
    if len(carriers) < 10:
        return dict(breadth=np.nan)
    # background of each carrier
    bgs = {}
    for s in carriers:
        v = int(s)
        while par[v] != -1 and (t[int(par[v])] - t[v]) < 1e-3:
            v = int(par[v])
        while par[v] != -1 and t[int(par[v])] <= tau:
            v = int(par[v])
        bgs[v] = bgs.get(v, 0) + 1
    n = len(carriers)
    broad = [b for b, c in bgs.items() if c >= 0.05 * n]
    dom = max(bgs, key=bgs.get)
    # background size = number of samples in that clade
    bsize = 0
    stack = [dom]
    while stack:
        u = stack.pop(); bsize += 1
        stack.extend(children.get(u, []))
    return dict(breadth=len(broad), dom_share=round(bgs[dom] / n, 3),
                dom_freq=round(bgs[dom] / max(bsize, 1), 3), carriers=n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bed')
    ap.add_argument('--pop', default='EUR')
    ap.add_argument('--breadth', action='store_true')
    ap.add_argument('--tau', type=float, default=200.0)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    bed = read_bed(a.bed)
    print(f'{len(bed)} loci from {a.bed}')
    rows = []
    scores = {}
    for r in bed.itertuples():
        if r.chrom not in scores:
            scores[r.chrom] = load_scores(r.chrom, a.pop)
        rec = dict(chrom=r.chrom, start=r.start, end=r.end, name=r.name,
                   segdup=is_segdup(r.chrom, (r.start + r.end) // 2))
        rec.update(score_region(r, scores[r.chrom]))
        if a.breadth and rec.get('n_variants') and np.isfinite(rec.get('best_pos', np.nan)):
            try:
                rec.update(breadth_at(r.chrom, rec['best_pos'], a.tau, a.pop))
            except Exception as e:
                rec['breadth'] = np.nan
        rows.append(rec)
    res = pd.DataFrame(rows)
    out = a.out or (OUT / 'loci_scan.csv')
    res.to_csv(out, index=False)
    print(f'wrote {out}')
    pd.set_option('display.width', 200)
    cols = [c for c in ('name', 'chrom', 'n_variants', 'k5', 'expected', 'enrichment',
                        'poisson_p', 'best_pct', 'best_age', 'segdup', 'breadth', 'dom_share') if c in res.columns]
    print(res[cols].to_string(index=False))


if __name__ == '__main__':
    main()
