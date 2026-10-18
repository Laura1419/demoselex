#!/usr/bin/env python3
"""1000G high-coverage phased SNPs -> ARG -> clade-age screen. Same analysis as before,
different substrate: 3,202 people (6,404 real haplotypes) instead of 503 imputed EUR.

Streams the phased VCF with cyvcf2, adds only common SNV sites to a transient tsinfer
SampleData (no big matrix in RAM, no retained VCF), infers, dumps trees, deletes the
intermediate. GRCh38 coordinates.

Stages: build <chrom> | date <chrom> | stats <chrom>
"""
import os, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/g1k'
OUT.mkdir(parents=True, exist_ok=True)
URL = ('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
       '1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/'
       '1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz')
MIN_MAF = 0.05
CHUNK = 10_000_000
HALF = 250_000


def vcf_path(chrom):
    return OUT / f'chr{chrom}.vcf.gz'


def download(chrom):
    """Resumable download (curl -C -); returns when the file is complete."""
    import subprocess
    p = vcf_path(chrom)
    done = p.with_suffix('.done')
    url = URL.format(chrom=chrom)
    if done.exists():
        print(f'chr{chrom}: download already complete', flush=True)
        return p
    t0 = time.time()
    for attempt in range(5):
        subprocess.run(['curl', '-sSL', '-C', '-', '-o', str(p), url], check=False)
        # verify against Content-Length
        head = subprocess.run(['curl', '-sI', url], capture_output=True, text=True).stdout
        want = None
        for line in head.splitlines():
            if line.lower().startswith('content-length:'):
                want = int(line.split(':')[1].strip())
        have = p.stat().st_size if p.exists() else 0
        if want and have == want:
            done.write_text('ok')
            print(f'chr{chrom}: downloaded {have/2**20:.0f} MB in {time.time()-t0:.0f}s', flush=True)
            return p
        print(f'chr{chrom}: partial {have}/{want} (attempt {attempt+1}), resuming', flush=True)
    raise RuntimeError(f'chr{chrom}: download did not complete')


def build(chrom):
    """Stream the VCF into a tsinfer SampleData, then infer and dump trees."""
    import cyvcf2, tsinfer
    trees = OUT / f'chr{chrom}.trees'
    if trees.exists():
        print(f'chr{chrom}: trees exist, skip'); return
    p = download(chrom)
    vcf = cyvcf2.VCF(str(p))
    sd_path = OUT / f'chr{chrom}.samples'
    if sd_path.exists():
        sd_path.unlink()
    seqlen = 260_000_000 if chrom == 1 else 250_000_000
    sd = tsinfer.SampleData(path=str(sd_path), sequence_length=seqlen)
    for _ in vcf.samples:
        sd.add_individual(ploidy=2, time=0.0)
    n_sites = n_seen = 0
    last_pos = -1
    t0 = time.time()
    for rec in vcf:
        n_seen += 1
        if rec.POS <= last_pos:
            continue                       # duplicated position (multiallelic split)
        if len(rec.REF) != 1 or any(len(a) != 1 for a in rec.ALT):
            continue                                   # SNVs only
        af = rec.INFO.get('AF')
        if isinstance(af, (list, tuple)):
            af = af[0] if len(af) else None
        if af is None:
            continue
        if min(af, 1 - af) < MIN_MAF:
            continue
        g = np.fromiter((a for gt in rec.genotypes for a in gt[:2]), dtype=np.int8,
                        count=2 * len(vcf.samples))
        sd.add_site(position=rec.POS, genotypes=g, alleles=[rec.REF] + list(rec.ALT))
        last_pos = rec.POS
        n_sites += 1
        if n_sites % 50000 == 0:
            print(f'  {n_sites} sites kept / {n_seen} seen ({time.time()-t0:.0f}s)', flush=True)
    sd.finalise()
    print(f'chr{chrom}: {n_sites} sites kept of {n_seen} records, '
          f'SampleData {sd_path.stat().st_size/2**20:.0f} MB', flush=True)
    sd.close()
    t1 = time.time()
    ts = tsinfer.infer(tsinfer.load(str(sd_path)), path_compression=True, progress_monitor=False)
    ts.dump(str(trees))
    print(f'chr{chrom}: trees {ts.num_trees} sites {ts.num_sites} samples {ts.num_samples} '
          f'in {time.time()-t1:.0f}s; trees size {trees.stat().st_size/2**20:.0f} MB', flush=True)
    sd_path.unlink(missing_ok=True)
    p.unlink(missing_ok=True)
    (OUT / f'chr{chrom}.done').unlink(missing_ok=True)   # do not keep the VCF


def date(chrom):
    import tskit, tsdate
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
        t0 = time.time()
        dated = tsdate.date(sub, mutation_rate=1.25e-8, population_size=10_000,
                            method='inside_outside', time_units='generations')
        dated.dump(str(dst))
        print(f'  chr{chrom} {start}-{end}: {sub.num_sites} sites, {time.time()-t0:.0f}s', flush=True)


def stats(chrom):
    import tskit, glob, pandas as pd
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
            rows.append((int(p), n, age, rho))
    pd.DataFrame(rows, columns=['pos', 'carriers', 'age', 'rho']).to_csv(
        OUT / f'chr{chrom}_stats.csv', index=False)
    print(f'chr{chrom}: {len(rows)} variants scored', flush=True)


def main():
    mode = sys.argv[1]
    chrom = int(sys.argv[2])
    if mode == 'build':
        build(chrom)
    elif mode == 'download':
        download(chrom)
    elif mode == 'date':
        date(chrom)
    elif mode == 'stats':
        stats(chrom)


if __name__ == '__main__':
    main()
