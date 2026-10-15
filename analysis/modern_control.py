#!/usr/bin/env python3
"""Test A: modern control (EUR1KG 503 samples, same imputed bundle).

Builds 1 Mb windows (same QC) for the modern samples, runs the validated ELAI
configuration (S=2, K=10, 100 steps, --pkk-site), and performs the tract test:
mean run-length of each individual's state along chr1 vs a shuffled-state null.

Writes data/derived/modern_control/.
"""
import json, re, subprocess, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
from load_akbari_geno import AkbariGeno, iter_pvar, parse_psam, paths  # noqa: E402

OUT = ROOT / 'data/derived/modern_control'
BINARY = '/home/erikg/elai-fork/src/elai'
S, K, STEPS, MG, SEED = 2, 10, 100, 10, 20261020
WINDOWS = [(1, 50_000_000), (1, 100_000_000), (1, 150_000_000),
           (1, 55_000_000), (1, 105_000_000), (1, 155_000_000),
           (1, 60_000_000), (1, 110_000_000), (1, 160_000_000),
           (1, 65_000_000), (1, 115_000_000), (1, 165_000_000),
           (1, 70_000_000), (1, 120_000_000), (1, 170_000_000),
           (1, 75_000_000), (1, 125_000_000), (1, 175_000_000),
           (1, 80_000_000), (1, 130_000_000), (1, 180_000_000),
           (1, 85_000_000), (1, 135_000_000), (1, 185_000_000),
           (1, 90_000_000), (1, 140_000_000), (1, 190_000_000),
           (1, 95_000_000), (1, 145_000_000), (1, 195_000_000)]


def marker_qc(h):
    valid = (h >= 0) & (h <= 2)
    calls = valid.mean(1)
    af = np.divide(np.where(valid, h, 0).sum(1), 2 * valid.sum(1),
                   out=np.zeros(len(h), dtype=float), where=valid.sum(1) > 0)
    maf = np.minimum(af, 1 - af)
    return (calls >= .95) & (maf >= .05)


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    ids = [str(x) for x in parse_psam(paths()[2])]
    kg = [i for i in ids if re.match(r'^(NA|HG)\d{5}$', i)]
    idx = {sid: k for k, sid in enumerate(ids)}
    sub = np.array([idx[s] for s in kg], dtype=np.uint32)
    json.dump(kg, open(OUT / 'sample_ids.json', 'w'))
    print(f'modern samples: {len(kg)}', flush=True)
    wins = {}
    with AkbariGeno(sample_subset=sub) as r, open(OUT / 'build.log', 'w') as lg:
        for c, posc in WINDOWS:
            lo, hi = posc - 500_000, posc + 500_000
            batch = []
            for i, (ch, p, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
                if ch == str(c) and lo <= p <= hi and len(ref) == 1 and len(alt) == 1 \
                        and ref in 'ACGT' and alt in 'ACGT' and ref != alt:
                    batch.append((i, p, vid, ref, alt))
            if not batch:
                continue
            vi = np.array([b[0] for b in batch], dtype=np.uint32)
            h = r.read_hardcalls(vi)
            mask = marker_qc(h)
            keep = [(b, h[j]) for j, b in enumerate(batch) if mask[j]]
            if len(keep) < 150:
                print(f'  c{c}_w{posc//10**6}: only {len(keep)} QC markers, skip', flush=True)
                continue
            if len(keep) > 800:
                step = len(keep) / 800
                keep = [keep[int(k * step)] for k in range(800)]
            tag = f'c{c}_w{posc // 10**6}'
            wd = OUT / tag
            wd.mkdir(exist_ok=True)
            center = posc
            with (wd / 'genotypes.txt').open('w') as g:
                for (i, p, vid, ref, alt), gen in keep:
                    codes = np.array([ref + ref, ref + alt, alt + alt, '??'])
                    gv = np.where(gen > 2, 3, gen)
                    g.write(vid + ' ' + ' '.join(codes[gv]) + '\n')
            with (wd / 'panel.inp').open('w') as f:
                f.write(f"{len(kg)}\n{len(keep)}\nIND " + ' '.join(kg) + '\n')
                with (wd / 'genotypes.txt').open() as src:
                    f.write(src.read())
            with (wd / 'panel.pos').open('w') as f:
                for (i, p, vid, ref, alt), gen in keep:
                    f.write(f'{vid} {p} {c}\n')
            wins[tag] = dict(chrom=c, center=center, markers=len(keep))
            print(f'  {tag}: {len(keep)} markers', flush=True)
    json.dump(wins, open(OUT / 'windows.json', 'w'))


def run():
    wins = json.load(open(OUT / 'windows.json'))
    n = len(json.load(open(OUT / 'sample_ids.json')))
    for tag in wins:
        wd = OUT / tag
        prefix = f'w_s{S}_k{K}_seed{SEED}'
        done = wd / 'output' / f'{prefix}.pkk-site.tsv'
        if done.exists():
            continue
        cmd = [BINARY, '-g', 'panel.inp', '-p', '1', '-pos', 'panel.pos', '-C', str(S), '-c', str(K),
               '-w', '0', '-s', str(STEPS), '-mg', str(MG), '-R', str(SEED), '-o', prefix,
               '--pkk-site', f"{wins[tag]['chrom']}:{wins[tag]['center']}", '--pk1']
        with (wd / f'{prefix}.launch.log').open('w') as log:
            p = subprocess.run(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
        ok = p.returncode == 1 and '-elai: finished' in (wd / f'{prefix}.launch.log').read_text()
        if ok:
            import pandas as pd
            d = pd.read_csv(done, sep='\t')
            kcols = K * (K + 1) // 2
            d = d[d.k1 <= d.k2].copy()
            d['u'] = d.k1.astype(str) + '|' + d.k2.astype(str)
            q = d.groupby(['individual_index', 'u']).q_unordered.sum().unstack(fill_value=0.0)
            Q = np.zeros((n, kcols))
            kk = [(a, b) for a in range(K) for b in range(a, K)]
            cmap = {f'{a}|{b}': s for s, (a, b) in enumerate(kk)}
            for c_ in q.columns:
                if c_ in cmap:
                    Q[:, cmap[c_]] = q[c_].to_numpy()
            np.save(OUT / tag / f'{prefix}.q.npy', Q.astype(np.float32))
        print(f'{tag}: {"ok" if ok else "FAIL"}', flush=True)


def tract_test():
    wins = sorted(json.load(open(OUT / 'windows.json')))
    import pandas as pd
    samples = json.load(open(OUT / 'sample_ids.json'))
    n = len(samples)
    M = np.zeros((n, len(wins)), dtype=np.int64)
    C = np.zeros((n, len(wins)), dtype=bool)
    for j, t in enumerate(wins):
        f = OUT / t / 'output' / f'w_s{S}_k{K}_seed{SEED}.pk1.txt'
        # pk1 gives per-haplotype lower-cluster marginals; use hard argmax of diploid sum
        d = pd.read_csv(f, sep=r'\s+', header=None).to_numpy() if f.exists() else None
        if d is None:
            continue
    print('pk1 not needed: use q.npy path instead', flush=True)
    # simpler: use --pkk-site q matrices (1000.. no, 503 x 55)
    for j, t in enumerate(wins):
        f = OUT / t / f'w_s{S}_k{K}_seed{SEED}.q.npy'
        if not f.exists():
            continue
        Q = np.load(f)
        M[:, j] = Q.argmax(1)
        C[:, j] = Q.max(1) > 0.9
    def runlen(x):
        if len(x) < 3:
            return np.nan
        return len(x) / (np.sum(x[1:] != x[:-1]) + 1)
    rng = np.random.default_rng(7)
    real, shuf = [], []
    for i in range(n):
        x = M[i][C[i]]
        if len(x) < 10:
            continue
        real.append(runlen(x))
        vals = []
        for _ in range(30):
            y = x.copy()
            rng.shuffle(y)
            vals.append(runlen(y))
        shuf.append(np.median(vals))
    real, shuf = np.array(real), np.array(shuf)
    ratio = real / shuf
    print(f'\nTRACT TEST — modern 1000G, {len(real)} individuals x {len(wins)} windows:')
    print(f'  mean real run length: {np.nanmean(real):.2f} windows | shuffled null: {np.nanmean(shuf):.2f}')
    print(f'  ratio real/shuffled: median {np.nanmedian(ratio):.2f} | >1.5: {100*np.nanmean(ratio>1.5):.0f}% of individuals')
    conf_frac = C.mean(1)
    print(f'  mean confident windows per individual: {conf_frac.mean():.1f} of {len(wins)}')
    json.dump(dict(n_individuals=int(len(real)),
                   mean_run_real=float(np.nanmean(real)), mean_run_shuf=float(np.nanmean(shuf)),
                   median_ratio=float(np.nanmedian(ratio)),
                   frac_ratio_gt1_5=float(np.nanmean(ratio > 1.5)),
                   mean_confident_windows=float(conf_frac.mean()),
                   n_windows=len(wins)),
              open(OUT / 'tract_test.json', 'w'), indent=1)


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if mode in ('all', 'build'):
        build()
    if mode in ('all', 'run'):
        run()
    if mode in ('all', 'test'):
        tract_test()


if __name__ == '__main__':
    main()
