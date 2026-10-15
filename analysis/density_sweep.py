#!/usr/bin/env python3
"""Density sweep with a label-consistent design (fixes the run-length confound).

One ELAI fit per (cohort, density) over a contiguous 4 Mb block; a ladder of 8
bridge sites at 500 kb spacing inside the block, each exported from the SAME
fitted model (same seed -> identical EM trajectory -> consistent state labels
across sites, verified by em.txt identity). Run-length of each individual's
state across the ladder is then a VALID tract statistic, unlike the previous
cross-independent-run test.

Cohorts: modern (503 1000G WGS samples, real high-density genotypes) and
ancient (1,000 dated, capture-imputed). Densities: 800/Mb, 2000/Mb, 5000/Mb
markers, uniform-spaced within the block after the standard QC.

Writes data/derived/density_sweep/.
"""
import json, subprocess, sys, zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, iter_pvar, paths, parse_psam  # noqa: E402

OUT = ROOT / 'data/derived/density_sweep'
BINARY = '/home/erikg/elai-fork/src/elai'
S, K, STEPS, MG = 2, 10, 100, 10
SEED = 20261020
SEED_B = 20261021
CHROM, LO, HI = 1, 50_000_000, 54_000_000
DENSITIES = {'d800': 800, 'd2000': 2000, 'd5000': 5000}   # markers per Mb
SITES = [LO + 250_000 + 500_000 * i for i in range(8)]    # 8 ladder sites
ANCIENT_SAMPLES = ROOT / 'data/derived/elai_amy_dense/work_long/samples.json'
MIN_MARKERS = 300


def marker_qc(h):
    valid = (h >= 0) & (h <= 2)
    calls = valid.mean(1)
    af = np.divide(np.where(valid, h, 0).sum(1), 2 * valid.sum(1),
                   out=np.zeros(len(h), dtype=float), where=valid.sum(1) > 0)
    maf = np.minimum(af, 1 - af)
    return (calls >= .95) & (maf >= .05)


def collect(cohort):
    """Return (indices array, sample ids) for the cohort."""
    if cohort == 'ancient':
        samples = json.load(open(ANCIENT_SAMPLES))
        return np.array([s['index'] for s in samples], dtype=np.uint32), [s['sample'] for s in samples]
    ids = [str(x) for x in parse_psam(paths()[2])]
    kg = [i for i in ids if re.match(r'^(NA|HG)\d{5}$', i)]
    idx = {sid: k for k, sid in enumerate(ids)}
    return np.array([idx[s] for s in kg], dtype=np.uint32), kg


def build(cohort):
    vi_all, names = collect(cohort)
    rows = []
    with AkbariGeno(sample_subset=vi_all) as r:
        batch = []
        for i, (ch, p, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
            if ch == str(CHROM) and LO <= p <= HI and len(ref) == 1 and len(alt) == 1 \
                    and ref in 'ACGT' and alt in 'ACGT' and ref != alt:
                batch.append((i, p, vid, ref, alt))
        # process in chunks to bound memory
        kept = []
        for k in range(0, len(batch), 4096):
            chunk = batch[k:k + 4096]
            vi = np.array([b[0] for b in chunk], dtype=np.uint32)
            h = r.read_hardcalls(vi)
            mask = marker_qc(h)
            kept.extend((b, h[j].copy()) for j, b in enumerate(chunk) if mask[j])
    print(f'{cohort}: {len(kept)} QC markers in block', flush=True)
    for dname, per_mb in DENSITIES.items():
        cap = per_mb * (HI - LO) // 1_000_000
        tag = f'{cohort}_{dname}'
        wd = OUT / tag
        wd.mkdir(parents=True, exist_ok=True)
        kept.sort(key=lambda b: b[0][1])
        if len(kept) > cap:
            step = len(kept) / cap
            kept_d = [kept[int(k * step)] for k in range(cap)]
        else:
            kept_d = kept
        if len(kept_d) < MIN_MARKERS:
            print(f'{tag}: too few markers ({len(kept_d)})'); continue
        with (wd / 'genotypes.txt').open('w') as g:
            for (i, p, vid, ref, alt), gen in kept_d:
                codes = np.array([ref + ref, ref + alt, alt + alt, '??'])
                gv = np.where(gen > 2, 3, gen)
                g.write(vid + ' ' + ' '.join(codes[gv]) + '\n')
        with (wd / 'panel.inp').open('w') as f:
            f.write(f"{len(names)}\n{len(kept_d)}\nIND " + ' '.join(names) + '\n')
            with (wd / 'genotypes.txt').open() as src:
                f.write(src.read())
        with (wd / 'panel.pos').open('w') as f:
            for (i, p, vid, ref, alt), gen in kept_d:
                f.write(f'{vid} {p} {CHROM}\n')
        json.dump(dict(markers=len(kept_d), sites=SITES), open(wd / 'meta.json', 'w'))
        print(f'  {tag}: {len(kept_d)} markers written', flush=True)


def run_one(tag, site, seed):
    wd = OUT / tag
    prefix = f'w_s{S}_k{K}_seed{seed}_s{site}'
    out = wd / 'output' / f'{prefix}.pkk-site.tsv'
    if out.exists():
        return True
    cmd = [BINARY, '-g', 'panel.inp', '-p', '1', '-pos', 'panel.pos', '-C', str(S), '-c', str(K),
           '-w', '0', '-s', str(STEPS), '-mg', str(MG), '-R', str(seed), '-o', prefix,
           '--pkk-site', f'{CHROM}:{site}']
    with (wd / f'{prefix}.launch.log').open('w') as log:
        p = subprocess.run(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT, timeout=14400)
    ok = p.returncode == 1 and '-elai: finished' in (wd / f'{prefix}.launch.log').read_text()
    return ok


def q_from_site(tag, site, seed, n):
    wd = OUT / tag
    prefix = f'w_s{S}_k{K}_seed{seed}_s{site}'
    d = pd.read_csv(wd / 'output' / f'{prefix}.pkk-site.tsv', sep='\t')
    assert abs(d.groupby('individual_index').p_ordered.sum().median() - 1) < 1e-5
    d = d[d.k1 <= d.k2].copy()
    d['u'] = d.k1.astype(str) + '|' + d.k2.astype(str)
    q = d.groupby(['individual_index', 'u']).q_unordered.sum().unstack(fill_value=0.0)
    kk = [(a, b) for a in range(K) for b in range(a, K)]
    cmap = {f'{a}|{b}': s for s, (a, b) in enumerate(kk)}
    Q = np.zeros((n, len(kk)))
    for c_ in q.columns:
        if c_ in cmap:
            Q[:, cmap[c_]] = q[c_].to_numpy()
    return Q


def run(tag, seed=SEED):
    meta = json.load(open(OUT / tag / 'meta.json'))
    n = len(open(OUT / tag / 'panel.inp').readline().split()) - 1
    n = int(open(OUT / tag / 'panel.inp').readline())
    for site in meta['sites']:
        ok = run_one(tag, site, seed)
        print(f'{tag} s{site}: {"ok" if ok else "FAIL"}', flush=True)


def test(cohort):
    rows = []
    n_map = {'ancient': 1000, 'modern': None}
    for dname in DENSITIES:
        tag = f'{cohort}_{dname}'
        wd = OUT / tag
        if not (wd / 'meta.json').exists():
            continue
        n = int(open(wd / 'panel.inp').readline())
        Qs = {}
        for site in SITES:
            f = wd / 'output' / f'w_s{S}_k{K}_seed{SEED}_s{site}.pkk-site.tsv'
            if f.exists():
                Qs[site] = q_from_site(tag, site, SEED, n)
        if len(Qs) < 4:
            continue
        sites = sorted(Qs)
        H = np.stack([Qs[s].argmax(1) for s in sites], 1)     # (n, sites)
        Cmax = np.stack([Qs[s].max(1) for s in sites], 1)
        # label-consistency check: same-seed fits must be identical up to site; verify via
        # em.txt equality across two site fits
        ems = [ (wd / 'output' / f'w_s{S}_k{K}_seed{SEED}_s{sites[0]}.em.txt').read_text(),
                (wd / 'output' / f'w_s{S}_k{K}_seed{SEED}_s{sites[-1]}.em.txt').read_text()]
        labels_consistent = ems[0] == ems[1]
        # run length per individual (confident sites only)
        def runlen(x):
            x = x[~np.isnan(x)]
            if len(x) < 3:
                return np.nan
            return len(x) / (np.sum(x[1:] != x[:-1]) + 1)
        rng = np.random.default_rng(11)
        real, shuf = [], []
        same_frac = []
        for i in range(n):
            x = H[i].astype(float)
            ok = Cmax[i] > 0.9
            x[~ok] = np.nan
            if np.isfinite(x).sum() < 4:
                continue
            real.append(runlen(x))
            same = (x[1:] == x[:-1]) & np.isfinite(x[1:]) & np.isfinite(x[:-1])
            if np.isfinite(x[:-1]).sum() >= 2:
                same_frac.append(np.nanmean(same))
            vals = []
            f = x[~np.isnan(x)].copy()
            for _ in range(50):
                rng.shuffle(f)
                vals.append(runlen(f))
            shuf.append(np.median(vals))
        real, shuf = np.array(real), np.array(shuf)
        rec = dict(cohort=cohort, density=dname, n_ind=len(real),
                   mean_run=float(np.nanmean(real)), mean_shuf=float(np.nanmean(shuf)),
                   ratio=float(np.nanmean(real) / np.nanmean(shuf)),
                   median_ratio=float(np.nanmedian(real / shuf)),
                   frac_same_adjacent=float(np.nanmean(same_frac)),
                   labels_consistent=bool(labels_consistent),
                   mean_conf=float(np.nanmean(Cmax > 0.9)))
        rows.append(rec)
        print(json.dumps(rec), flush=True)
    return rows


import re  # noqa: E402


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    OUT.mkdir(parents=True, exist_ok=True)
    if mode in ('all', 'build'):
        for c in ('modern', 'ancient'):
            build(c)
    if mode in ('all', 'run'):
        import concurrent.futures as cf
        jobs = []
        for c in ('modern', 'ancient'):
            for d in DENSITIES:
                tag = f'{c}_{d}'
                if (OUT / tag / 'meta.json').exists():
                    jobs.append(tag)
        with cf.ThreadPoolExecutor(max_workers=8) as ex:
            list(ex.map(lambda t: run(t), jobs))
    if mode in ('all', 'test'):
        allr = []
        for c in ('modern', 'ancient'):
            allr += test(c)
        pd.DataFrame(allr).to_csv(OUT / 'density_sweep_results.csv', index=False)
        print('WROTE', OUT / 'density_sweep_results.csv')


if __name__ == '__main__':
    main()
