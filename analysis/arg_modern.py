#!/usr/bin/env python3
"""Modern-only ARG (contemporaneous samples) + frequency trajectory.

Samples are all time-0 (1000G EUR), so the lineage-fraction trajectory is not
confounded by temporally structured sampling, and tsdate's inside-outside method
(which requires contemporaneous samples) applies. The inferred rise window is then
cross-checked against the ancient allele frequencies across time.

Usage: arg_modern.py <chrom> <lo> <hi> [label:pos,...]
"""
import json, re, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, paths, parse_psam  # noqa: E402

OUT = ROOT / 'data/derived/argmodern'
OUT.mkdir(parents=True, exist_ok=True)
MAP = ROOT / 'data/derived/argmerge/map_all.json'


def modern_indices():
    ids = [str(x) for x in parse_psam(paths()[2])]
    idx, names = [], []
    for i, s in enumerate(ids):
        if re.match(r'^(NA|HG)\d{5}$', s):
            idx.append(i); names.append(s)
    return np.array(idx, dtype=np.uint32), names


def build(chrom, lo, hi, max_sites=40000, min_maf=0.05):
    mp = [m for m in json.load(open(MAP))[str(chrom)] if lo <= m['pos'] <= hi]
    mp.sort(key=lambda m: m['pos'])
    if len(mp) > max_sites:
        st = len(mp) / max_sites
        mp = [mp[int(i * st)] for i in range(max_sites)]
    idx, names = modern_indices()
    reader = AkbariGeno(sample_subset=idx)
    var = np.array([m['akb_var'] for m in mp], dtype=np.uint32)
    rng = np.random.default_rng(20261020)
    G = np.empty((len(mp), len(idx)), dtype=np.int8)
    for k in range(0, len(mp), 2048):
        h = reader.read_hardcalls(var[k:k + 2048]).astype(np.int8)
        g = np.full_like(h, -1)
        g[h == 0] = 0
        g[h == 2] = 1
        het = h == 1                      # unphased heterozygote: resolve at random
        g[het] = rng.integers(0, 2, size=int(het.sum())).astype(np.int8)
        G[k:k + len(h)] = g
    called = (G >= 0).mean(1)
    n_alt = (G == 1).sum(1); n_ref = (G == 0).sum(1)
    af = n_alt / np.maximum(n_alt + n_ref, 1)
    keep = (called >= 0.9) & (np.minimum(af, 1 - af) >= min_maf) & (n_alt > 0) & (n_ref > 0)
    G = G[keep]; pos = np.array([m['pos'] for m in mp])[keep]
    np.savez_compressed(OUT / f'chr{chrom}_{lo}_{hi}.npz', G=G, pos=pos, names=np.array(names))
    print(f'chr{chrom}:{lo}-{hi}: {G.shape[0]} sites x {G.shape[1]} moderns', flush=True)


def infer(chrom, lo, hi):
    import tsinfer, time as _t
    d = np.load(OUT / f'chr{chrom}_{lo}_{hi}.npz', allow_pickle=True)
    G, pos = d['G'], d['pos']
    rel = (pos - lo).astype(np.int32)
    p = OUT / f'chr{chrom}_{lo}_{hi}.samples'
    if p.exists():
        p.unlink()
    sd = tsinfer.SampleData(path=str(p), sequence_length=int(rel.max()) + 1_000_000)
    for _ in range(G.shape[1]):
        sd.add_individual(ploidy=1, time=0.0)
    for i in range(G.shape[0]):
        g = G[i]
        sd.add_site(position=int(rel[i]), genotypes=g, alleles=['A', 'T'])
    sd.finalise(); sd.close()
    sd = tsinfer.load(str(p))
    t0 = _t.time()
    ts = tsinfer.infer(sd, path_compression=True, progress_monitor=False)
    ts.dump(str(OUT / f'chr{chrom}_{lo}_{hi}.trees'))
    print(f'inferred: trees {ts.num_trees} sites {ts.num_sites} {_t.time()-t0:.0f}s', flush=True)
    p.unlink(missing_ok=True)


def date(chrom, lo, hi, Ne=10_000):
    import tskit, tsdate, time as _t
    ts = tskit.load(str(OUT / f'chr{chrom}_{lo}_{hi}.trees')).simplify()
    t0 = _t.time()
    dated = tsdate.date(ts, mutation_rate=1.25e-8, population_size=Ne,
                        method='inside_outside', time_units='generations')
    dated.dump(str(OUT / f'chr{chrom}_{lo}_{hi}.dated.trees'))
    print(f'dated {_t.time()-t0:.0f}s', flush=True)


def main():
    mode, chrom, lo, hi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    if mode == 'build':
        build(chrom, lo, hi)
    elif mode == 'infer':
        infer(chrom, lo, hi)
    elif mode == 'date':
        date(chrom, lo, hi)
    elif mode == 'traj':
        import tskit
        import importlib
        at = importlib.import_module('arg_traj')
        ts = tskit.load(str(OUT / f'chr{chrom}_{lo}_{hi}.dated.trees'))
        print(f'ARG sites {ts.num_sites} samples {ts.num_samples} trees {ts.num_trees}')
        targets = [(s.split(':')[0], int(s.split(':')[1])) for s in sys.argv[5].split(',')]
        for lab, p in targets:
            r = at.analyse(ts, p - lo, lab)
            if r is None:
                print(f'{lab}: no mutation'); continue
            print(f"{lab}: freq={r['freq']:.3f} origin={r['origin_gen']:.0f} max_rho={r['max_rho']:.3f} "
                  f"T_rise={r.get('T_rise', float('nan')):.0f} width={r.get('width', float('nan')):.0f}")
            c = r['curve']
            for i in range(0, len(c), max(1, len(c) // 12)):
                print(f'    {c[i][0]:9.1f} gen : {c[i][1]:.3f}')


if __name__ == '__main__':
    main()
