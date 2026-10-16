#!/usr/bin/env python3
"""Within-tree permutation null for ARG selection candidates.

For a variant with n carriers, the null is: draw random sets of n samples FROM THE SAME
LOCAL TREE and compute the same statistic. This preserves demography, recombination
context and local Ne (they are baked into the genealogy) and tests exactly one thing:
is this carrier set unusually clumped / recent for this tree? No simulation, no
demographic model.

Statistics per set S of samples in a local tree:
  mrca_time  = time of the most recent common ancestor of S
  span_len   = total length of the minimal subtree connecting S (clade compactness)

Usage:
  arg_permnull.py <chrom> <N|targets.csv> [K]
    N            score the top-N variants by our z on that chromosome
    targets.csv  a csv with columns chrom,pos
"""
import sys, json
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/arggenome'
PERM = 200


def stats_for_set(par, t, nodes):
    """MRCA time and minimal-subtree length for a set of nodes."""
    counts = {}
    for u in nodes:
        v = u
        while v != -1:
            counts[v] = counts.get(v, 0) + 1
            v = par[v]
    n = len(nodes)
    mrca = None
    for v, c in counts.items():
        if c == n and (mrca is None or t[v] < t[mrca]):
            mrca = v
    if mrca is None:
        return np.nan, np.nan
    span = 0.0
    for v in counts:
        if v == mrca:
            continue
        p = par[v]
        if p != -1:
            span += t[p] - t[v]
    return float(t[mrca]), float(span)


def run(chrom, targets, K=PERM, seed=20261020):
    rng = np.random.default_rng(seed)
    rows = []
    by_chunk = {}
    for r in targets.itertuples():
        by_chunk.setdefault(int(r.chunk), []).append((int(r.pos), float(r.z)))
    for chunk, items in sorted(by_chunk.items()):
        f = OUT / f'chr{chrom}_{chunk}_{chunk + 10_000_000}.dated.trees'
        if not f.exists():
            print(f'  missing {f.name}, skipping {len(items)} targets'); continue
        ts = tskit.load(str(f))
        t = ts.tables.nodes.time
        G = ts.genotype_matrix()
        spos = np.array([s.position for s in ts.sites()])
        samples = ts.samples()
        print(f'  chunk {chunk}: {len(items)} targets, {ts.num_sites} sites', flush=True)
        for pos, z_existing in items:
            rel = pos - chunk
            j = int(np.argmin(np.abs(spos - rel)))
            if abs(spos[j] - rel) > 1:
                continue
            carriers = samples[G[j] == 1]
            n = len(carriers)
            if n < 3:
                continue
            tree = ts.at(spos[j])
            par = tree.parent_array
            obs_m, obs_s = stats_for_set(par, t, carriers)
            nm, ns = [], []
            for _ in range(K):
                rs = rng.choice(samples, size=n, replace=False)
                m, s_ = stats_for_set(par, t, rs)
                nm.append(m); ns.append(s_)
            nm = np.array(nm); ns = np.array(ns)
            def pval(null, obs, side='lower'):
                if side == 'lower':
                    return float((np.sum(null <= obs) + 1) / (K + 1))
                return float((np.sum(null >= obs) + 1) / (K + 1))
            rows.append(dict(
                chrom=chrom, pos=int(pos), n=n,
                mrca=obs_m, mrca_null_med=float(np.median(nm)), p_mrca=pval(nm, obs_m),
                z_mrca=float((obs_m - nm.mean()) / nm.std()) if nm.std() > 0 else np.nan,
                span=obs_s, span_null_med=float(np.median(ns)), p_span=pval(ns, obs_s),
                z_span=float((obs_s - ns.mean()) / ns.std()) if ns.std() > 0 else np.nan,
                z_heuristic=z_existing))
    return pd.DataFrame(rows)


def main():
    chrom = int(sys.argv[1])
    spec = sys.argv[2]
    K = int(sys.argv[3]) if len(sys.argv) > 3 else PERM
    z = pd.read_csv(OUT / 'genomewide_z.csv')
    if spec.isdigit():
        sub = z[z.chrom == chrom].nsmallest(int(spec), 'z')
    else:
        sub = pd.read_csv(spec)
    df = run(chrom, sub, K=K)
    out = OUT / f'permnull_chr{chrom}_{spec if not spec.isdigit() else "top" + spec}.csv'
    df.to_csv(out, index=False)
    print(f'wrote {out} ({len(df)} variants)')
    if len(df):
        print('\ncombined p (max of the two) < 0.05: %d | < 0.01: %d' %
              (int(((df[["p_mrca", "p_span"]].max(1)) < .05).sum()),
               int(((df[["p_mrca", "p_span"]].max(1)) < .01).sum())))
        print('\nstrongest by permutation null (p_mrca):')
        print(df.nsmallest(10, 'p_mrca')[['pos', 'n', 'mrca', 'mrca_null_med', 'p_mrca', 'p_span', 'z_heuristic']].to_string(index=False))


if __name__ == '__main__':
    main()
