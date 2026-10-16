#!/usr/bin/env python3
"""Merged-substrate ARG: AADR 1240K ancients (real calls) + Akbari-bundle 1000G moderns,
matched site-by-site on rsID. The only combination that carries the canonical transition
variants (LCT rs4988235, SLC45A2, SLC24A5) AND present-day clock anchors.

Modes:
  map <chrom>    write the rsID intersection table for a chromosome
  build <chrom>  assemble the (sites x samples) haploid matrix -> npz
  infer <chrom>  tsinfer -> trees (radiocarbon ages attached)
  date <chrom>   tsdate -> dated trees
  test <chrom>   derived-allele age vs frequency + named target variants

Writes data/derived/argmerge/.
"""
import json, re, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_aadr import GenotypeMatrix, parse_ind, parse_anno, iter_snp  # noqa: E402
from load_akbari_geno import AkbariGeno, iter_pvar, paths, parse_psam  # noqa: E402

OUT = ROOT / 'data/derived/argmerge'
OUT.mkdir(parents=True, exist_ok=True)
AADR = ROOT / 'data/raw/aadr-v66.p1'
K1240 = AADR / 'v66.p1_1240K.aadr.patch.PUB'
ANCIENT_SAMPLES = ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'
YEARS_PER_GEN = 29.0
EUR_MODERN = {'Spanish', 'French', 'Basque', 'Sardinian', 'Icelandic', 'Norwegian', 'English'}


def rsid_key(x):
    x = str(x)
    m = re.search(r'(rs\d+)', x)
    return m.group(1) if m else None


def map_all():
    """Match 1240K and Akbari sites by (chrom, position); record allele orientation.

    Akbari variant IDs are chrom_pos_ref_alt (no rsIDs), so position is the join key.
    A site is usable when REF/ALT agree; when they are swapped the modern dosage is
    flipped at build time (flip=True).
    """
    by_chrom = {}
    for j, r in enumerate(iter_snp(str(K1240) + '.snp')):
        c = r[1]
        if not c.isdigit() or not 1 <= int(c) <= 22:
            continue
        by_chrom.setdefault(c, {})[int(r[3])] = (j, r[4], r[5], r[0])
    out, n_seen = {}, 0
    akb = {}
    for i, (c, p, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
        if c not in by_chrom:
            continue
        hit = by_chrom[c].get(p)
        if hit is None:
            continue
        j, ref0, alt0, rs = hit
        if (ref, alt) == (ref0, alt0):
            flip = False
        elif (ref, alt) == (alt0, ref0):
            flip = True
        else:
            continue
        out.setdefault(c, []).append(dict(aadr_row=j, pos=p, akb_var=i, flip=flip, rs=rs))
        n_seen += 1
    json.dump(out, open(OUT / 'map_all.json', 'w'))
    print('shared sites per chromosome:',
          {k: len(v) for k, v in sorted(out.items(), key=lambda kv: int(kv[0]))})
    return out


def collect_samples():
    ancients = json.load(open(ANCIENT_SAMPLES))
    ids = [str(x) for x in parse_psam(paths()[2])]
    aadr = list(parse_ind(str(K1240) + '.ind'))
    bybase = {}
    for i, ind in enumerate(aadr):
        bybase.setdefault(re.split(r'[._]', ind.id)[0], []).append(i)
    A, atimes, anames = [], [], []
    for a in ancients:
        b = re.split(r'[._]', a['sample'])[0]
        cand = bybase.get(b, [])
        if len(cand) != 1:
            continue
        A.append(cand[0]); atimes.append(a['BP'] / YEARS_PER_GEN); anames.append(a['sample'])
    # moderns: European groups with no radiocarbon date
    anno = parse_anno(str(AADR / 'v66.p1_1240K.aadr.PUB.anno'))

    def dated(i):
        r = anno.get(i.id)
        return r is not None and r.date_bp is not None
    # moderns come from the Akbari psam directly: AADR 1240K has no modern samples.
    # The bundle's EUR1KG set is 1000G WGS (real calls), IDs NAxxxxx / HGxxxxx.
    M, mnames = [], []
    for i, sid in enumerate(ids):
        if re.match(r'^(NA|HG)\d{5}$', sid):
            M.append(i); mnames.append(sid)
    print(f'ancients {len(A)} | moderns {len(M)}')
    return np.array(A, dtype=np.int64), np.array(atimes), anames, M, mnames


def build(chrom, min_calls=0.5, min_maf=0.03, max_sites=100000):
    mp = json.load(open(OUT / 'map_all.json'))[str(chrom)]
    a_idx, a_times, a_names, m_idx, m_names = collect_samples()
    n_a, n_m = len(a_idx), len(m_idx)
    rows = np.array([m['aadr_row'] for m in mp], dtype=np.int64)
    apos = np.array([m['pos'] for m in mp], dtype=np.int64)
    akb_var = np.array([m['akb_var'] for m in mp], dtype=np.uint32)
    # ancients (haploid from AADR): 0->0, 2->1, else missing
    gm = GenotypeMatrix(str(K1240) + '.geno')
    A = np.full((len(mp), n_a), -1, dtype=np.int8)
    for k, row in enumerate(rows):
        codes = gm.individual_genotypes(int(row), np.arange(gm.n_snp, dtype=np.int64)) \
            if False else None
    # faster: decode only the needed sites per individual
    for k in range(n_a):
        codes = gm.individual_genotypes(int(a_idx[k]), rows)
        g = np.full(len(mp), -1, dtype=np.int8)
        g[codes == 0] = 0
        g[codes == 2] = 1
        A[:, k] = g
        if (k + 1) % 250 == 0:
            print(f'  ancients {k+1}/{n_a}', flush=True)
    # moderns from the Akbari pgen (diploid -> haploid: 0->0, 2->1, het/missing -> -1)
    reader = AkbariGeno(sample_subset=np.array(m_idx, dtype=np.uint32))
    Mp = np.empty((len(mp), n_m), dtype=np.int8)
    for k in range(0, len(mp), 2048):
        chunk = akb_var[k:k + 2048]
        h = reader.read_hardcalls(chunk).astype(np.int8)
        g = np.full_like(h, -1)
        g[h == 0] = 0
        g[h == 2] = 1
        Mp[k:k + len(chunk)] = g
    G = np.concatenate([A, Mp], axis=1)
    called = (G >= 0).mean(1)
    n_alt = (G == 1).sum(1); n_ref = (G == 0).sum(1)
    af = np.divide(n_alt, np.maximum(n_alt + n_ref, 1))
    keep = (called >= min_calls) & (np.minimum(af, 1 - af) >= min_maf) & (n_alt > 0) & (n_ref > 0)
    G, pos = G[keep], apos[keep]
    if len(pos) > max_sites:
        st = len(pos) / max_sites
        js = [int(i * st) for i in range(max_sites)]
        G, pos = G[js], pos[js]
    times = np.concatenate([a_times, np.zeros(n_m)])
    names = a_names + m_names
    np.savez_compressed(OUT / f'chr{chrom}_geno.npz', G=G, pos=pos,
                        times=times.astype(np.float32), names=np.array(names))
    print(f'chr{chrom}: {G.shape[0]} sites x {G.shape[1]} samples', flush=True)


def infer(chrom):
    import tsinfer
    d = np.load(OUT / f'chr{chrom}_geno.npz', allow_pickle=True)
    G, pos, times = d['G'], d['pos'], d['times']
    p = OUT / f'chr{chrom}.samples'
    if p.exists():
        p.unlink()
    sd = tsinfer.SampleData(path=str(p), sequence_length=int(pos.max()) + 1_000_000)
    for t in times:
        sd.add_individual(ploidy=1, time=float(t))
    for i in range(G.shape[0]):
        g = G[i]
        if (g == 1).sum() == 0 or (g == 0).sum() == 0:
            continue
        sd.add_site(position=int(pos[i]), genotypes=g, alleles=['A', 'T'])
    sd.finalise(); sd.close()
    sd = tsinfer.load(str(p))
    import time as _t
    t0 = _t.time()
    ts = tsinfer.infer(sd, path_compression=True, progress_monitor=False)
    tbl = ts.dump_tables()
    t = tbl.nodes.time.copy() + float(times.max()) + 1.0
    smp = ts.samples()
    t[smp] = times
    tbl.nodes.time = t
    tbl.sort()
    tbl.tree_sequence().dump(str(OUT / f'chr{chrom}.trees'))
    print(f'chr{chrom}: trees {ts.num_trees} sites {ts.num_sites} {_t.time()-t0:.0f}s', flush=True)
    p.unlink(missing_ok=True)


def date(chrom):
    import tskit, tsdate, time as _t
    ts = tskit.load(str(OUT / f'chr{chrom}.trees'))
    t0 = _t.time()
    ts = ts.simplify()          # tsdate requires no unary nodes
    dated = tsdate.date(ts, mutation_rate=1.25e-8, time_units='generations')
    dated.dump(str(OUT / f'chr{chrom}.dated.trees'))
    print(f'chr{chrom}: dated {_t.time()-t0:.0f}s', flush=True)


def test(chrom, targets=None):
    import tskit, pandas as pd
    ts = tskit.load(str(OUT / f'chr{chrom}.dated.trees'))
    time = ts.tables.nodes.time
    mut = ts.tables.mutations
    pos = np.array([s.position for s in ts.sites()])
    Gm = ts.genotype_matrix()
    freq = Gm.sum(1) / Gm.shape[1]
    age = np.full(ts.num_sites, np.nan)
    for s, nd in zip(mut.site, mut.node):
        t = time[nd]
        if not np.isfinite(age[s]) or t < age[s]:
            age[s] = t
    df = pd.DataFrame(dict(pos=pos, freq=freq, age=age))
    d = df[(df.freq > 0) & (df.age > 0)].copy()
    d['bin'] = pd.qcut(d.freq, 10, duplicates='drop')
    curve = d.groupby('bin', observed=True).age.agg(['median', 'std', 'count'])
    print(curve.to_string())
    if targets:
        for lab, p in targets:
            m = df[np.abs(df.pos - p) <= 1]
            if not len(m):
                print(f'{lab} @{p}: ABSENT from ARG sites'); continue
            r = m.iloc[0]
            sub = d[d.bin == r.bin] if False else None
            med, sd = curve.loc[pd.cut([r.freq], curve.index.categories)[0], ['median', 'std']] \
                if hasattr(curve.index, 'categories') else (np.nan, np.nan)
            z = (r.age - med) / sd if sd and sd > 0 else np.nan
            print(f'{lab} @{int(r.pos)} freq={r.freq:.3f} age={r.age:.1f} bin_med={med:.1f} z={z:.2f}')


def main():
    mode = sys.argv[1]
    chrom = int(sys.argv[2])
    if mode == 'map':
        map_all()
    elif mode == 'build':
        build(chrom)
    elif mode == 'infer':
        infer(chrom)
    elif mode == 'date':
        date(chrom)
    elif mode == 'test':
        t = None
        if len(sys.argv) > 3:
            t = [(s.split('@')[0], int(s.split('@')[1])) for s in sys.argv[3].split(',')]
        test(chrom, t)


if __name__ == '__main__':
    main()
