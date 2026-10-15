#!/usr/bin/env python3
"""Dense-substrate ARG (Akbari imputed bundle) for a target region.

Builds a tsinfer ARG over a genomic window from the dense 'AllVariants' pgen
(ancients + modern EUR), dates it with tsdate using radiocarbon ages, and runs the
derived-allele-age-vs-frequency test at that window. This is the substrate with
enough sites to resolve allele ages, and it contains the canonical transition
variants (LCT rs4988235 etc.) that HO lacks.

Usage: arg_dense.py <chrom> <start> <end> [label]
"""
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import tskit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, iter_pvar, paths, parse_psam  # noqa: E402

OUT = ROOT / 'data/derived/argdense'
OUT.mkdir(parents=True, exist_ok=True)
YEARS_PER_GEN = 29.0
ANCIENT_SAMPLES = ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'


def collect_samples():
    ancients = json.load(open(ANCIENT_SAMPLES))
    ids = [str(x) for x in parse_psam(paths()[2])]
    idx = {s: i for i, s in enumerate(ids)}
    subs, times, names, kinds = [], [], [], []
    for a in ancients:
        subs.append(a['index']); times.append(a['BP'] / YEARS_PER_GEN)
        names.append(a['sample']); kinds.append('ancient')
    import re
    modern = [i for i in ids if re.match(r'^(NA|HG)\d{5}$', i)]
    for m in modern:
        subs.append(idx[m]); times.append(0.0); names.append(m); kinds.append('modern')
    order = np.argsort(np.array(subs, dtype=np.uint32))
    subs = np.array(subs, dtype=np.uint32)[order]
    times = np.array(times, dtype=float)[order]
    names = [names[i] for i in order]
    kinds = [kinds[i] for i in order]
    return subs, times, names, kinds


def build_region(chrom, start, end, min_maf=0.05, max_sites=20000):
    subs, times, names, kinds = collect_samples()
    rows = []
    for i, (c, p, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
        if c == str(chrom) and start <= p <= end and len(ref) == 1 and len(alt) == 1 \
                and ref in 'ACGT' and alt in 'ACGT' and ref != alt:
            rows.append((i, p, vid, ref, alt))
    print(f'region sites: {len(rows)}', flush=True)
    reader = AkbariGeno(sample_subset=subs)
    G = np.zeros((len(rows), len(subs)), dtype=np.int8)
    for k in range(0, len(rows), 4096):
        chunk = rows[k:k + 4096]
        h = reader.read_hardcalls(np.array([b[0] for b in chunk], dtype=np.uint32)).astype(np.int8)
        h[h > 2] = -1
        G[k:k + len(chunk)] = h
    called = (G >= 0).mean(1)
    af = np.divide(np.where(G >= 0, G, 0).sum(1), 2 * np.maximum((G >= 0).sum(1), 1))
    maf = np.minimum(af, 1 - af)
    keep = (called >= 0.5) & (maf >= min_maf)
    G, rows = G[keep], [r for r, k in zip(rows, keep) if k]
    if len(rows) > max_sites:
        step = len(rows) / max_sites
        idx = [int(i * step) for i in range(max_sites)]
        G, rows = G[idx], [rows[i] for i in idx]
    print(f'kept {len(rows)} sites (called>=0.5, maf>={min_maf})', flush=True)
    np.savez_compressed(OUT / f'c{chrom}_{start}_{end}_geno.npz', G=G,
                        pos=np.array([r[1] for r in rows]),
                        vid=np.array([r[2] for r in rows]),
                        times=times, names=np.array(names), kinds=np.array(kinds))
    return len(rows)


def infer(chrom, start, end):
    import tsinfer
    d = np.load(OUT / f'c{chrom}_{start}_{end}_geno.npz', allow_pickle=True)
    G, pos, times = d['G'], d['pos'], d['times']
    path = OUT / f'c{chrom}_{start}_{end}.samples'
    if path.exists():
        path.unlink()
    sd = tsinfer.SampleData(path=str(path), sequence_length=int(pos.max()) + 1_000_000)
    for t in times:
        sd.add_individual(ploidy=1, time=float(t))
    for i in range(G.shape[0]):
        g = G[i]
        if (g == 1).sum() == 0 or (g == 0).sum() == 0:
            continue
        sd.add_site(position=int(pos[i]), genotypes=g, alleles=['A', 'T'], inference=None)
    sd.finalise()
    sd.close()
    print('SampleData written', flush=True)
    sd = tsinfer.load(str(path))
    ts = tsinfer.infer(sd, path_compression=True, progress_monitor=False)
    # attach ages: shift internals above oldest sample, keep the topology valid
    tbl = ts.dump_tables()
    t = tbl.nodes.time.copy() + float(times.max()) + 1.0
    smp = ts.samples()
    t[smp] = times
    tbl.nodes.time = t
    tbl.sort()
    ts = tbl.tree_sequence()
    ts.dump(str(OUT / f'c{chrom}_{start}_{end}.trees'))
    print('trees', ts.num_trees, 'sites', ts.num_sites, flush=True)


def date(chrom, start, end):
    import tsdate
    ts = tskit.load(str(OUT / f'c{chrom}_{start}_{end}.trees'))
    dated = tsdate.date(ts, mutation_rate=1.25e-8, time_units='generations')
    dated.dump(str(OUT / f'c{chrom}_{start}_{end}.dated.trees'))
    print('dated', flush=True)


def test(chrom, start, end, variants=None):
    ts = tskit.load(str(OUT / f'c{chrom}_{start}_{end}.dated.trees'))
    time = ts.tables.nodes.time
    mut = ts.tables.mutations
    pos = np.array([s.position for s in ts.sites()])
    G = ts.genotype_matrix()
    freq = G.sum(1) / G.shape[1]
    age = np.full(ts.num_sites, np.nan)
    for s, nd in zip(mut.site, mut.node):
        t = time[nd]
        if not np.isfinite(age[s]) or t < age[s]:
            age[s] = t
    df = pd.DataFrame(dict(pos=pos, freq=freq, age=age))
    d = df[(df.freq > 0) & (df.age > 0)].copy()
    d['bin'] = pd.qcut(d.freq, 8, duplicates='drop')
    curve = d.groupby('bin', observed=True).age.agg(['median', 'std', 'count'])
    print(curve.to_string())
    if variants:
        for vid, p in variants:
            m = df[np.abs(df.pos - p) < 1]
            if len(m):
                r = m.iloc[0]
                b = pd.cut([r.freq], curve.index.categories if hasattr(curve.index, 'categories') else None)
                print(f'{vid} pos={int(r.pos)} freq={r.freq:.3f} age={r.age:.1f}')
            else:
                print(f'{vid} pos={p}: not in ARG sites')


def main():
    mode, chrom, start, end = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    if mode == 'build':
        build_region(chrom, start, end)
    elif mode == 'infer':
        infer(chrom, start, end)
    elif mode == 'date':
        date(chrom, start, end)
    elif mode == 'test':
        test(chrom, start, end)


if __name__ == '__main__':
    main()
