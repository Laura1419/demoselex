#!/usr/bin/env python3
"""ARG screen, rebuilt properly: REAL PHASE + DENSE sites.

Two things were wrong before:
  (1) genotypes were read unphased and heterozygotes were coin-flipped into pseudo-
      haplotypes, fabricating haplotype sharing that does not exist;
  (2) the site universe was the 1240K x Akbari intersection (~1 site / 2.7 kb), discarding
      ~98% of available variants.
The bundle carries a phase bit (verified: 2,514 / 2,515 calls) and 9.7M high-quality
variants, so both are fixed here: 503 EUR moderns -> 1,006 real haplotypes, sites from
HighQualityVariants.tsv.

Stages
  index                 map HighQualityVariants -> pgen variant indices (cached)
  build <chrom>         phased haplotype matrix for the chromosome -> npz
  infer <chrom>         tsinfer over 1,006 haplotypes -> .trees
  date <chrom>          tsdate inside-outside in 10 Mb chunks
  stats <chrom>         carrier-clade age + rho + compartment -> csv

Writes data/derived/argphased/.
"""
import gzip, json, re, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, paths, parse_psam  # noqa: E402

OUT = ROOT / 'data/derived/argphased'
OUT.mkdir(parents=True, exist_ok=True)
HQV = ROOT / 'data/raw/akbari/tier2/HighQualityVariants.tsv'
SEGDUP = ROOT / 'data/raw/annotation/segdup/genomicSuperDups.txt.gz'
CHUNK = 10_000_000
MIN_MAF = 0.05
HALF = 250_000


def modern_haps():
    """Indices of the 503 EUR modern samples in the psam."""
    ids = [str(x) for x in parse_psam(paths()[2])]
    return np.array([i for i, s in enumerate(ids) if re.match(r'^(NA|HG)\d{5}$', s)], dtype=np.uint32)


def hqv_keys():
    keys = set()
    with open(HQV) as f:
        for line in f:
            if line.startswith('#'):
                continue
            p = line.split('\t')
            if len(p) >= 2:
                keys.add(f'{p[0]}_{p[1]}')
    print(f'HighQualityVariants: {len(keys):,} sites', flush=True)
    return keys


def index():
    """(chrom, pos) -> pgen variant index, restricted to the high-quality set."""
    keys = hqv_keys()
    out = {str(c): ([], []) for c in range(1, 23)}
    from load_akbari_geno import iter_pvar
    for i, (c, p, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
        if c not in out:
            continue
        if f'{c}_{p}' in keys and len(ref) == 1 and len(alt) == 1 and ref in 'ACGT' and alt in 'ACGT':
            out[c][0].append(i); out[c][1].append(p)
        if i % 5_000_000 == 0:
            print(f'  scanned {i:,}', flush=True)
    np.savez_compressed(OUT / 'hqv_index.npz',
                        **{f'idx_{c}': np.array(v[0], dtype=np.int64) for c, v in out.items()},
                        **{f'pos_{c}': np.array(v[1], dtype=np.int64) for c, v in out.items()})
    print('sites per chromosome:', {c: len(v[0]) for c, v in out.items() if v[0]}, flush=True)


def load_index(chrom):
    z = np.load(OUT / 'hqv_index.npz')
    return z[f'idx_{chrom}'], z[f'pos_{chrom}']


def build(chrom):
    if (OUT / f'chr{chrom}_hap.npz').exists():
        print(f'chr{chrom}: exists, skip'); return
    idx, pos = load_index(chrom)
    samp = modern_haps()
    n_h = 2 * len(samp)
    reader = AkbariGeno(sample_subset=samp)
    pr = reader._reader
    H = np.full((len(idx), n_h), -1, dtype=np.int8)
    al = np.empty(n_h, dtype=np.int32)
    ph = np.empty(len(samp), dtype=np.int8)
    for k, v in enumerate(idx):
        pr.read_alleles_and_phasepresent(int(v), al, ph)
        H[k] = np.where(al <= 1, al, -1)          # 0=ref, 1=alt, else missing
        if (k + 1) % 50000 == 0:
            print(f'  {k+1}/{len(idx)}', flush=True)
    called = (H >= 0).mean(1)
    af = np.where(called > 0, (H == 1).sum(1) / np.maximum((H >= 0).sum(1), 1), 0.0)
    keep = (called >= 0.9) & (np.minimum(af, 1 - af) >= MIN_MAF)
    H, pos = H[keep], pos[keep]
    order = np.argsort(pos, kind='stable')
    H, pos = H[order], pos[order]
    uniq = np.concatenate([[True], np.diff(pos) > 0])     # drop repeated positions
    H, pos = H[uniq], pos[uniq]
    np.savez_compressed(OUT / f'chr{chrom}_hap.npz', H=H, pos=pos)
    print(f'chr{chrom}: {H.shape[0]} sites x {H.shape[1]} haplotypes (kept from {len(idx)})', flush=True)


def infer(chrom):
    import tsinfer, time as _t
    if (OUT / f'chr{chrom}.trees').exists():
        print(f'chr{chrom}: trees exist, skip'); return
    d = np.load(OUT / f'chr{chrom}_hap.npz')
    H, pos = d['H'], d['pos']
    p = OUT / f'chr{chrom}.samples'
    if p.exists():
        p.unlink()
    sd = tsinfer.SampleData(path=str(p), sequence_length=int(pos.max()) + 1_000_000)
    for _ in range(H.shape[1]):
        sd.add_individual(ploidy=1, time=0.0)
    for i in range(H.shape[0]):
        g = H[i]
        if (g == 1).sum() == 0 or (g == 0).sum() == 0:
            continue
        sd.add_site(position=int(pos[i]), genotypes=g, alleles=['A', 'T'])
    sd.finalise(); sd.close()
    t0 = _t.time()
    ts = tsinfer.infer(tsinfer.load(str(p)), path_compression=True, progress_monitor=False)
    ts.dump(str(OUT / f'chr{chrom}.trees'))
    print(f'chr{chrom}: trees {ts.num_trees} sites {ts.num_sites} in {_t.time()-t0:.0f}s', flush=True)
    p.unlink(missing_ok=True)


def date(chrom):
    import tskit, tsdate, time as _t
    ts = tskit.load(str(OUT / f'chr{chrom}.trees'))
    for start in range(0, int(ts.sequence_length), CHUNK):
        end = min(start + CHUNK, int(ts.sequence_length))
        dst = OUT / f'chr{chrom}_{start}_{end}.dated.trees'
        if dst.exists():
            continue
        if not any(start <= s.position < end for s in ts.sites()):
            continue
        sub = ts.keep_intervals([[start, end]]).trim().simplify()
        if sub.num_sites == 0:
            continue
        t0 = _t.time()
        dated = tsdate.date(sub, mutation_rate=1.25e-8, population_size=10_000,
                            method='inside_outside', time_units='generations')
        dated.dump(str(dst))
        print(f'  chr{chrom} {start}-{end}: {sub.num_sites} sites, {_t.time()-t0:.0f}s', flush=True)


MASK = None


def segdup_mask():
    global MASK
    if MASK is not None:
        return MASK
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
    out = {}
    for c, lst in iv.items():
        lst.sort(); merged = []
        for a, b in lst:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        out[c] = (np.array([m[0] for m in merged]), np.array([m[1] for m in merged]))
    MASK = out
    return out


def is_segdup(chrom, pos):
    m = segdup_mask().get(chrom)
    if m is None:
        return False
    s, e = m
    i = np.searchsorted(s, pos, side='right') - 1
    return bool(i >= 0 and pos <= e[i])


def stats(chrom):
    import tskit, glob
    rows = []
    for f in sorted(glob.glob(str(OUT / f'chr{chrom}_*.dated.trees'))):
        start = int(f.split('_')[-2])
        ts = tskit.load(f)
        t = ts.tables.nodes.time
        G = ts.genotype_matrix()
        spos = np.array([s.position for s in ts.sites()]) + start
        bps = np.array(list(ts.breakpoints())) + start
        node = {}
        for m in ts.mutations():
            k = int(ts.site(m.site).position) + start
            node.setdefault(k, m.node)
        for j, p in enumerate(spos):
            n = int(G[j].sum())
            if n < 3 or int(p) not in node:
                continue
            age = float(t[node[int(p)]])
            lo = np.searchsorted(bps, p - HALF); hi = np.searchsorted(bps, p + HALF)
            rho = (hi - lo) / (2 * HALF / 1000.0)
            rows.append((int(p), n, age, rho, 1 if is_segdup(chrom, int(p)) else 0))
    import pandas as pd
    pd.DataFrame(rows, columns=['pos', 'carriers', 'age', 'rho', 'segdup']).to_csv(
        OUT / f'chr{chrom}_stats.csv', index=False)
    print(f'chr{chrom}: {len(rows)} variants scored', flush=True)


def main():
    mode = sys.argv[1]
    if mode == 'index':
        index()
    elif mode == 'build':
        build(int(sys.argv[2]))
    elif mode == 'infer':
        infer(int(sys.argv[2]))
    elif mode == 'date':
        date(int(sys.argv[2]))
    elif mode == 'stats':
        stats(int(sys.argv[2]))


if __name__ == '__main__':
    main()
