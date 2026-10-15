#!/usr/bin/env python3
"""ARG pipeline on the raw AADR data (tsinfer + tsdate).

Substrate: HO (Human Origins) 1240K-adjacent array, transversion-heavy -> robust to
ancient DNA damage. Samples: our 1,000 radiocarbon-dated cohort (matched into AADR by
base identifier). Genotypes are the RAW (unimputed) AADR calls.

Modes
  map      match cohort ids into the AADR .ind, report coverage, save index map
  samples  build tsinfer SampleData for one chromosome (haploid sites) -> .samples
  infer    tsinfer.infer -> .trees
  date     tsdate.date with radiocarbon ages (generations) -> .dated.trees
  sanity   tree-sequence and time-distribution checks
  locus    carrier-clade / TMRCA test at a locus vs matched controls

Writes data/derived/arg/.
"""
import json, re, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_aadr import GenotypeMatrix, parse_ind, iter_snp, parse_anno  # noqa: E402

OUT = ROOT / 'data/derived/arg'
AADR = ROOT / 'data/raw/aadr-v66.p1'
HO = AADR / 'v66.p1_HO.aadr.patch.PUB'
K1240 = AADR / 'v66.p1_1240K.aadr.patch.PUB'
COHORT = ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'
YEARS_PER_GEN = 29.0


def base(sid):
    return re.split(r'[._]', sid)[0]


def load_cohort():
    return json.load(open(COHORT))


def build_map(dataset='HO'):
    """Match cohort sample ids to AADR individuals by base id."""
    cohort = load_cohort()
    prefix = HO if dataset in ('HO', 'mixed') else K1240
    inds = list(parse_ind(str(prefix) + '.ind'))
    by_base = {}
    for i, ind in enumerate(inds):
        by_base.setdefault(base(ind.id), []).append((i, ind.id))
    matched, ambiguous, missing = [], [], []
    for s in cohort:
        b = base(s['sample'])
        cand = by_base.get(b, [])
        if len(cand) == 1:
            matched.append(dict(cohort_id=s['sample'], aadr_index=cand[0][0], aadr_id=cand[0][1],
                                bp=s['BP'], country=s['country']))
        elif len(cand) == 0:
            missing.append(s['sample'])
        else:
            # prefer exact id match, else the one with matching .AG/.SG style
            exact = [c for c in cand if c[1] == s['sample']]
            if len(exact) == 1:
                matched.append(dict(cohort_id=s['sample'], aadr_index=exact[0][0], aadr_id=exact[0][1],
                                    bp=s['BP'], country=s['country']))
            else:
                ambiguous.append((s['sample'], [c[1] for c in cand]))
    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(dict(dataset=dataset, matched=matched, n_matched=len(matched),
                   n_ambiguous=len(ambiguous), n_missing=len(missing),
                   ambiguous=ambiguous[:40], missing=missing[:40]),
              open(OUT / f'map_{dataset}.json', 'w'), indent=1)
    print(json.dumps(dict(dataset=dataset, n_cohort=len(cohort), matched=len(matched),
                          ambiguous=len(ambiguous), missing=len(missing))))
    return matched


EUR_MODERN = {'Spanish', 'French', 'Basque', 'Sardinian', 'Icelandic', 'Norwegian', 'English'}


def build_map_mixed(base_dataset='HO', out_dataset='mixed'):
    """Ancients (from the validated map) + modern European samples at time 0.

    base_dataset selects the AADR bundle (HO = transversions only, 1240K = all sites);
    out_dataset names the emitted map/dataset."""
    prefix = HO if base_dataset == 'HO' else K1240
    ancients = json.load(open(OUT / f'map_{base_dataset}.json'))['matched']
    for a in ancients:
        a['cohort'] = 'ancient'
    inds = list(parse_ind(str(prefix) + '.ind'))
    anno = parse_anno(str(AADR / 'v66.p1_1240K.aadr.PUB.anno'))

    def dated(i):
        r = anno.get(i.id)
        return r is not None and r.date_bp is not None

    moderns = []
    for i, ind in enumerate(inds):
        if dated(ind) or ind.group not in EUR_MODERN:
            continue
        moderns.append(dict(cohort_id=ind.id, aadr_index=i, aadr_id=ind.id, bp=0.0,
                            country=ind.group, cohort='modern'))
    matched = ancients + moderns
    json.dump(dict(dataset=out_dataset, matched=matched, n_matched=len(matched),
                   n_ancient=len(ancients), n_modern=len(moderns)),
              open(OUT / f'map_{out_dataset}.json', 'w'), indent=1)
    print(json.dumps(dict(n_ancient=len(ancients), n_modern=len(moderns), total=len(matched))))


def chrom_sites(dataset, chrom):
    prefix = HO if dataset in ('HO', 'mixed') else K1240
    rows = []
    for j, r in enumerate(iter_snp(str(prefix) + '.snp')):
        if r[1] == str(chrom):
            rows.append((j, int(r[3]), r[0], r[4], r[5]))
    return rows


def build_samples(dataset='HO', chrom=1, min_calls=0.0):
    mapping = json.load(open(OUT / f'map_{dataset}.json'))['matched']
    rows = chrom_sites(dataset, chrom)
    n_sites = len(rows)
    n_ind = len(mapping)
    print(f'{dataset} chr{chrom}: {n_sites} sites x {n_ind} samples', flush=True)
    G = np.zeros((n_sites, n_ind), dtype=np.int8)   # 0=REF, 1=ALT, -1=missing
    site_idx = np.array([r[0] for r in rows], dtype=np.int64)
    prefix = HO if dataset in ('HO', 'mixed') else K1240
    reader = GenotypeMatrix(str(prefix) + ".geno")
    for k, m in enumerate(mapping):
        codes = reader.individual_genotypes(m['aadr_index'], site_idx)  # ALT dosage 0/1/2, 3 missing
        # haploid call: 0 -> ref(0), 2 -> alt(1), 1 (het from 2 reads) & 3 -> missing
        g = np.full(n_sites, -1, dtype=np.int8)
        g[codes == 0] = 0
        g[codes == 2] = 1
        G[:, k] = g
        if (k + 1) % 200 == 0:
            print(f'  decoded {k+1}/{n_ind}', flush=True)
    keep = (np.mean(G >= 0, axis=1) >= min_calls)
    # tsinfer needs at least one ALT and one REF per site
    n_alt = (G == 1).sum(1); n_ref = (G == 0).sum(1)
    keep &= (n_alt > 0) & (n_ref > 0)
    G = G[keep]; rows = [r for r, k in zip(rows, keep) if k]
    print(f'  kept {len(rows)} polymorphic sites', flush=True)
    np.savez_compressed(OUT / f'{dataset}_chr{chrom}_geno.npz',
                        G=G, pos=np.array([r[1] for r in rows]),
                        ids=np.array([m['cohort_id'] for m in mapping]),
                        bp=np.array([m['bp'] for m in mapping]))
    json.dump(dict(dataset=dataset, chrom=chrom, n_sites=len(rows), n_samples=n_ind,
                   mean_call=float(np.mean(json.load(open(OUT / f'map_{dataset}.json'))['matched'] and G >= 0))),
              open(OUT / f'{dataset}_chr{chrom}_sites.json', 'w'))
    return len(rows)


def build_tsinfer_samples(dataset='HO', chrom=1):
    import tsinfer
    d = np.load(OUT / f'{dataset}_chr{chrom}_geno.npz', allow_pickle=True)
    G, pos, ids, bp = d['G'], d['pos'], d['ids'], d['bp']
    chrom_len = int(pos.max()) + 1_000_000
    path = OUT / f'{dataset}_chr{chrom}.samples'
    if path.exists(): path.unlink()
    sd = tsinfer.SampleData(path=str(path), sequence_length=chrom_len)
    # individuals first (one haploid sample each), with radiocarbon ages
    for sid, y in zip(ids, bp):
        sd.add_individual(ploidy=1, time=float(y) / YEARS_PER_GEN)
    # then sites; genotypes columns are sample indices in order of addition
    for i in range(G.shape[0]):
        g = G[i]
        if (g == 1).sum() == 0 or (g == 0).sum() == 0:
            continue
        sd.add_site(position=int(pos[i]), genotypes=g, alleles=['A', 'T'])
    n_sites, n_samples = sd.num_sites, sd.num_samples
    sd.finalise()
    sd.close()
    print('wrote', path, 'sites', n_sites, 'samples', n_samples, flush=True)
    return path


def infer(dataset='HO', chrom=1, path_compression=False):
    import tsinfer
    import time as _t
    sd = tsinfer.load(str(OUT / f'{dataset}_chr{chrom}.samples'))
    t0 = _t.time()
    ts = tsinfer.infer(sd, path_compression=True, progress_monitor=False)
    # attach radiocarbon ages to the sample nodes. tsinfer leaves sample times at 0
    # (present day); our samples are ancient, so shift every internal node above the
    # oldest sample to keep the tree valid (parent strictly older than child). tsdate
    # overwrites all non-sample times, so the shift does not bias the dates.
    import numpy as _np
    bp = _np.array(json.load(open(OUT / f'map_{dataset}.json'))['matched'] and
                   [m['bp'] for m in json.load(open(OUT / f'map_{dataset}.json'))['matched']], dtype=float)
    sample_times = bp / YEARS_PER_GEN
    tbl = ts.dump_tables()
    times = tbl.nodes.time.copy()
    smp = ts.samples()
    assert len(smp) == len(sample_times), (len(smp), len(sample_times))
    offset = float(sample_times.max()) + 1.0
    times = times + offset
    times[smp] = sample_times
    tbl.nodes.time = times
    tbl.sort()
    ts = tbl.tree_sequence()
    print('sample times set: min %.1f max %.1f gen (offset %.1f)' %
          (ts.tables.nodes.time[ts.samples()].min(), ts.tables.nodes.time[ts.samples()].max(), offset),
          flush=True)
    out = OUT / f'{dataset}_chr{chrom}.trees'
    ts.dump(str(out))
    print('wrote', out, 'trees', ts.num_trees, 'edges', ts.num_edges,
          f'{_t.time()-t0:.0f}s', flush=True)


def date(dataset='HO', chrom=1, Ne=10_000):
    import tskit, tsdate, time as _t
    ts = tskit.load(str(OUT / f'{dataset}_chr{chrom}.trees'))
    ts = ts.simplify()
    t0 = _t.time()
    dated = tsdate.date(ts, mutation_rate=1.25e-8, time_units='generations')
    out = OUT / f'{dataset}_chr{chrom}.dated.trees'
    dated.dump(str(out))
    print('wrote', out, f'{_t.time()-t0:.0f}s', flush=True)


def sanity(dataset='HO', chrom=1):
    import tskit
    f = OUT / f'{dataset}_chr{chrom}.dated.trees'
    if not f.exists():
        f = OUT / f'{dataset}_chr{chrom}.trees'
    ts = tskit.load(str(f))
    ti = ts.tables.nodes.time
    samples = ts.samples()
    print(json.dumps(dict(file=f.name, trees=int(ts.num_trees), sites=int(ts.num_sites),
                          samples=int(ts.num_samples),
                          min_time=float(ti[samples].min()), max_time=float(ti.max()),
                          mean_tmrca=float(ti[samples].mean()))), flush=True)


def main():
    mode = sys.argv[1]
    ds = sys.argv[2] if len(sys.argv) > 2 else 'HO'
    ch = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    if mode == 'map':
        build_map(ds)
    elif mode == 'mapmix':
        build_map_mixed('HO', 'mixed')
    elif mode == 'map':
        build_map(ds)
    elif mode == 'mapmix1240':
        build_map('1240K')
        build_map_mixed('1240K', 'mixed1240')
    elif mode == 'samples':
        build_samples(ds, ch)
    elif mode == 'tsinfer':
        build_tsinfer_samples(ds, ch)
    elif mode == 'infer':
        infer(ds, ch)
    elif mode == 'date':
        date(ds, ch)
    elif mode == 'sanity':
        sanity(ds, ch)


if __name__ == '__main__':
    main()
