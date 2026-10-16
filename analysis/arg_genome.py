#!/usr/bin/env python3
"""Genome-wide ARG screen (moderns, merged substrate) — empirical z-score, no simulation.

Design (agreed): per-variant carrier-clade age from the dated ARG, standardised by an
empirical null computed within strata (frequency decile x local-recombination quartile),
ranked genome-wide. External calibration = known loci must land where they belong.

Stages (each independently parallelisable):
  build <chrom>    merged 1240K x Akbari sites for the chromosome, 503 EUR moderns,
                   unphased het resolved at random -> npz
  infer <chrom>    tsinfer -> .trees (samples all time 0)
  date <chrom>     tsdate inside-outside in 10 Mb chunks -> .dated.trees per chunk
  stats <chrom>    per-variant (pos, freq, clade age, local rho) -> csv
  pool             genome-wide strata, MAD z, rank, known-locus calibration

Writes data/derived/arggenome/.
"""
import json, re, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, paths, parse_psam  # noqa: E402

OUT = ROOT / 'data/derived/arggenome'
OUT.mkdir(parents=True, exist_ok=True)
MAP = ROOT / 'data/derived/argmerge/map_all.json'
CHUNK = 10_000_000
MIN_CALL, MIN_MAF = 0.9, 0.05
NE = 10_000


def modern_idx():
    ids = [str(x) for x in parse_psam(paths()[2])]
    return np.array([i for i, s in enumerate(ids) if re.match(r'^(NA|HG)\d{5}$', s)], dtype=np.uint32)


def build(chrom):
    if (OUT / f'chr{chrom}.npz').exists():
        print(f'chr{chrom}: npz exists, skip'); return
    mp = sorted([m for m in json.load(open(MAP))[str(chrom)]], key=lambda m: m['pos'])
    idx = modern_idx()
    reader = AkbariGeno(sample_subset=idx)
    var = np.array([m['akb_var'] for m in mp], dtype=np.uint32)
    rng = np.random.default_rng(20261020 + chrom)
    G = np.empty((len(mp), len(idx)), dtype=np.int8)
    for k in range(0, len(mp), 4096):
        h = reader.read_hardcalls(var[k:k + 4096]).astype(np.int8)
        g = np.full_like(h, -1)
        g[h == 0] = 0
        g[h == 2] = 1
        het = h == 1
        g[het] = rng.integers(0, 2, size=int(het.sum())).astype(np.int8)
        G[k:k + len(h)] = g
    called = (G >= 0).mean(1)
    n_alt = (G == 1).sum(1); n_ref = (G == 0).sum(1)
    af = n_alt / np.maximum(n_alt + n_ref, 1)
    keep = (called >= MIN_CALL) & (np.minimum(af, 1 - af) >= MIN_MAF) & (n_alt > 0) & (n_ref > 0)
    G = G[keep]
    pos = np.array([m['pos'] for m in mp])[keep]
    np.savez_compressed(OUT / f'chr{chrom}.npz', G=G, pos=pos)
    print(f'chr{chrom}: {G.shape[0]} sites x {G.shape[1]} moderns', flush=True)


def infer(chrom):
    import tsinfer, time as _t
    if (OUT / f'chr{chrom}.trees').exists():
        print(f'chr{chrom}: trees exist, skip'); return
    d = np.load(OUT / f'chr{chrom}.npz', allow_pickle=True)
    G, pos = d['G'], d['pos']
    p = OUT / f'chr{chrom}.samples'
    if p.exists():
        p.unlink()
    sd = tsinfer.SampleData(path=str(p), sequence_length=int(pos.max()) + 1_000_000)
    for _ in range(G.shape[1]):
        sd.add_individual(ploidy=1, time=0.0)
    for i in range(G.shape[0]):
        g = G[i]
        if (g == 1).sum() == 0 or (g == 0).sum() == 0:
            continue
        sd.add_site(position=int(pos[i]), genotypes=g, alleles=['A', 'T'])
    sd.finalise(); sd.close()
    sd = tsinfer.load(str(p))
    t0 = _t.time()
    ts = tsinfer.infer(sd, path_compression=True, progress_monitor=False)
    ts.dump(str(OUT / f'chr{chrom}.trees'))
    print(f'chr{chrom}: trees {ts.num_trees} sites {ts.num_sites} {_t.time()-t0:.0f}s', flush=True)
    p.unlink(missing_ok=True)


def date(chrom):
    import tskit, tsdate, time as _t
    src = OUT / f'chr{chrom}.trees'
    if not src.exists():
        print(f'chr{chrom}: no trees'); return
    ts = tskit.load(str(src))
    lo, hi = 0, int(ts.sequence_length)
    n = 0
    for start in range(lo, hi, CHUNK):
        end = min(start + CHUNK, hi)
        # only chunks that actually contain data
        if not any(start <= s.position < end for s in ts.sites()):
            continue
        dst = OUT / f'chr{chrom}_{start}_{end}.dated.trees'
        if dst.exists():
            n += 1; continue
        sub = ts.keep_intervals([[start, end]]).trim().simplify()
        if sub.num_sites == 0:
            continue
        t0 = _t.time()
        dated = tsdate.date(sub, mutation_rate=1.25e-8, population_size=NE,
                            method='inside_outside', time_units='generations')
        dated.dump(str(dst))
        n += 1
        print(f'  chr{chrom} {start}-{end}: sites {sub.num_sites} dated {_t.time()-t0:.0f}s', flush=True)
    print(f'chr{chrom}: {n} chunks dated', flush=True)


def stats(chrom):
    import tskit, glob
    if (OUT / f'chr{chrom}_stats.csv').exists():
        print(f'chr{chrom}: stats exist, skip'); return
    rows = []
    for f in sorted(glob.glob(str(OUT / f'chr{chrom}_*.dated.trees'))):
        start = int(f.split('_')[-2]); end = int(f.split('_')[-1].split('.')[0])
        ts = tskit.load(f)
        t = ts.tables.nodes.time
        mut = ts.tables.mutations
        pos = np.array([s.position for s in ts.sites()]) + start   # undo trim
        G = ts.genotype_matrix()
        freq = G.sum(1) / G.shape[1]
        age = np.full(ts.num_sites, np.nan)
        for s, nd in zip(mut.site, mut.node):
            v = t[nd]
            if not np.isfinite(age[s]) or v < age[s]:
                age[s] = v
        rho = ts.num_trees / max(end - start, 1)          # local recombination proxy
        for p, fr, ag in zip(pos, freq, age):
            if fr > 0 and np.isfinite(ag):
                rows.append((int(p), float(fr), float(ag), rho, int(start)))
    import pandas as pd
    df = pd.DataFrame(rows, columns=['pos', 'freq', 'age', 'rho', 'chunk'])
    df.to_csv(OUT / f'chr{chrom}_stats.csv', index=False)
    print(f'chr{chrom}: {len(df)} variants scored', flush=True)


def pool():
    import pandas as pd
    files = sorted(OUT.glob('chr*_stats.csv'))
    frames = []
    for f in files:
        chrom = int(f.name.split('_')[0].replace('chr', ''))
        t = pd.read_csv(f); t['chrom'] = chrom; frames.append(t)
    df = pd.concat(frames, ignore_index=True)
    print(f'pooled: {len(df)} variants from {len(files)} chromosomes')
    df['fbin'] = pd.qcut(df.freq, 10, duplicates='drop')
    df['rbin'] = pd.qcut(df.rho.rank(method='first'), 4, labels=False, duplicates='drop')
    g = df.groupby(['fbin', 'rbin'], observed=True).age
    med = g.transform('median'); mad = g.transform(lambda x: (x - x.median()).abs().median())
    df['z'] = (df.age - med) / (1.4826 * mad.replace(0, np.nan))
    df.to_csv(OUT / 'genomewide_z.csv', index=False)
    print('\ntop 25 by z (young clade for its frequency x recombination context):')
    cols = ['pos', 'freq', 'age', 'rho', 'z']
    print(df.nsmallest(25, 'z')[cols].to_string(index=False))
    known = {'LCT': 136608646, 'SLC45A2': 33951693, 'SLC24A5': 48426484}
    print('\nknown-locus calibration (rank among all variants, 1 = most extreme):')
    df['rank'] = df.z.rank(method='min')
    for lab, p in known.items():
        m = df[np.abs(df.pos - p) <= 1]
        if len(m):
            r = m.iloc[0]
            print(f'  {lab:9s} pos={int(r.pos)} freq={r.freq:.3f} age={r.age:.0f} z={r.z:.2f} '
                  f'rank={int(r["rank"])}/{len(df)} percentile={100*(df.z<r.z).mean():.2f}%')
        else:
            print(f'  {lab:9s} not in scored set')


def main():
    mode = sys.argv[1]
    if mode == 'build':
        build(int(sys.argv[2]))
    elif mode == 'infer':
        infer(int(sys.argv[2]))
    elif mode == 'date':
        date(int(sys.argv[2]))
    elif mode == 'stats':
        stats(int(sys.argv[2]))
    elif mode == 'pool':
        pool()


if __name__ == '__main__':
    main()
