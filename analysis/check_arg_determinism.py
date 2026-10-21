#!/usr/bin/env python3
"""Determinism check for the ARG build: infer the same window twice, compare exactly.

Needs the 1000G chr22 phased VCF (only ~425 MB). Builds a 2 Mb window through the identical
code path twice, runs tsinfer on each, and compares the tree sequences element by element
(edges, nodes, mutations, sites) - and, if they match, confirms tsdate output matches too.
"""
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import tskit

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
OUT = ROOT / 'data/derived/detcheck'
OUT.mkdir(parents=True, exist_ok=True)
VCF = OUT / 'chr22.vcf.gz'
URL = ('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/'
       '1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/'
       '1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz')
LO, HI = 10_000_000, 12_000_000


def ensure_vcf():
    if VCF.exists() and VCF.stat().st_size > 400_000_000:
        return
    print('downloading chr22 VCF (~425 MB)', flush=True)
    subprocess.run(['curl', '-sSL', '-C', '-', '-o', str(VCF), URL], check=True)
    print('downloaded', VCF.stat().st_size // 2**20, 'MB', flush=True)


def build_and_infer(tag):
    """Same code path as analysis/g1k_phased.py build(), restricted to a 2 Mb window."""
    import cyvcf2, tsinfer
    p = OUT / f'{tag}.samples'
    if p.exists():
        p.unlink()
    sd = tsinfer.SampleData(path=str(p), sequence_length=HI + 1_000_000)
    vcf = cyvcf2.VCF(str(VCF))
    for _ in vcf.samples:
        sd.add_individual(ploidy=2, time=0.0)
    last = -1
    kept = 0
    for rec in vcf:
        if not (LO <= rec.POS < HI):
            continue
        if rec.POS <= last:
            continue
        if len(rec.REF) != 1 or any(len(a) != 1 for a in rec.ALT):
            continue
        af = rec.INFO.get('AF')
        if isinstance(af, (list, tuple)):
            af = af[0] if len(af) else None
        if af is None or min(af, 1 - af) < 0.05:
            continue
        g = np.fromiter((a for gt in rec.genotypes for a in gt[:2]), dtype=np.int8,
                        count=2 * len(vcf.samples))
        sd.add_site(position=rec.POS, genotypes=g, alleles=[rec.REF] + list(rec.ALT))
        last = rec.POS
        kept += 1
    sd.finalise(); sd.close()
    print(f'{tag}: {kept} sites', flush=True)
    ts = tsinfer.infer(tsinfer.load(str(p)), path_compression=True, progress_monitor=False)
    p.unlink(missing_ok=True)
    return ts


def compare(a, b, label):
    out = {}
    out['num_nodes'] = (a.num_nodes, b.num_nodes)
    out['num_edges'] = (a.num_edges, b.num_edges)
    out['num_trees'] = (a.num_trees, b.num_trees)
    ok = True
    for name in ('edges', 'nodes', 'mutations', 'sites'):
        ta, tb = getattr(a.tables, name), getattr(b.tables, name)
        ok &= bool(ta == tb)
        out[name] = 'identical' if ta == tb else 'DIFFER'
        ok = ok and out[name] == 'identical'
    out['max_time_diff'] = float(np.max(np.abs(a.tables.nodes.time - b.tables.nodes.time)))
    out['identical'] = bool(ok)
    print(label, json.dumps(out), flush=True)
    return out


def main():
    ensure_vcf()
    a = build_and_infer('runA')
    b = build_and_infer('runB')
    r1 = compare(a, b, 'tsinfer A vs B (same input, same code):')
    if r1['identical']:
        import tsdate
        da = tsdate.date(a.simplify(), mutation_rate=1.25e-8, population_size=10_000,
                         method='inside_outside', time_units='generations')
        db = tsdate.date(b.simplify(), mutation_rate=1.25e-8, population_size=10_000,
                         method='inside_outside', time_units='generations')
        compare(da, db, 'tsdate A vs B:')
    json.dump(r1, open(OUT / 'determinism.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
