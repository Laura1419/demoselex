"""Site-ladder tract-continuity probe (S=2, K=10, 1 Mb window, 2 seeds).

Physical-locality test: export the local-state posterior at a ladder of sites
spanning the AMY region and check whether the posterior is piecewise-constant -
stable within ~tract-length stretches, switching between them - and whether the
same boundaries recur across independently seeded fits. Reads the marker list
from the existing 1 Mb build; writes only under data/derived/elai_ws/.

Sites: 10 positions at 300 kb spacing centred on the AMY virtual site
(chr1:104,121,000), i.e. chr1:102,721,000 .. 105,521,000. All inside the 1 Mb
build's marker span (chr1:103.62-104.62 Mb) must be replaced by bracketing
markers via the patch's missing-emission bridge, so ladder positions outside
the build are legal: the patch bridges to the nearest bracketing markers.
"""
import json, subprocess, sys
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

BINARY = '/home/erikg/elai-fork/src/elai'
BUILD = Path('data/derived/elai_ws/w1000kb')
OUTDIR = Path('data/derived/elai_ws/site_ladder')
SITE_BASE = 104121000
SITES = [103671000, 103821000, 104121000, 104421000, 104571000]
SEEDS = [20261020, 20261021]
K, S, STEPS, MG = 10, 2, 100, 10


def marker_span():
    positions = [m[2] for m in json.load(open(BUILD / '..', 'window_scaling.report.json'))] if False else None
    pos = np.loadtxt(BUILD / 'panel.pos', usecols=1, dtype=np.int64)
    return int(pos.min()), int(pos.max()), len(pos)


def run_site(site, seed):
    tag = f's{site // 1000}k'
    dst = OUTDIR / tag
    dst.mkdir(parents=True, exist_ok=True)
    # link the existing 1 Mb build inputs (read-only reuse)
    for f in ('panel.inp', 'panel.pos', 'genotypes.txt'):
        t = dst / f
        if not t.exists():
            t.symlink_to((BUILD / f).resolve())
    prefix = f'w_s{S}_k{K}_seed{seed}_site{site}'
    cmd = [BINARY, '-g', 'panel.inp', '-p', '1', '-pos', 'panel.pos', '-C', str(S), '-c', str(K),
           '-w', '0', '-s', str(STEPS), '-mg', str(MG), '-R', str(seed), '-o', prefix,
           '--pkk-site', f'1:{site}']
    with (dst / f'{prefix}.launch.log').open('w') as log:
        proc = subprocess.run(cmd, cwd=dst, stdout=log, stderr=subprocess.STDOUT, timeout=3600)
    ok = proc.returncode == 1 and '-elai: finished' in (dst / f'{prefix}.launch.log').read_text()
    return tag, prefix, ok


def diag(tag, prefix):
    d = pd.read_csv(OUTDIR / tag / 'output' / f'{prefix}.pkk-site.tsv', sep='\t')
    assert abs(d.groupby('individual_index').p_ordered.sum().median() - 1) < 1e-5
    mx = d.groupby('individual_index').p_ordered.max()
    marg = d.pivot_table(index='individual_index', columns='k1', values='p_ordered', aggfunc='sum')
    return mx.median(), (mx > .5).mean(), marg.idxmax(axis=1).to_numpy()


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    lo, hi, nmk = marker_span()
    print(f'1 Mb build marker span: chr1:{lo}-{hi} ({nmk} markers); ladder {SITES[0]}..{SITES[-1]}', flush=True)
    runs = []
    for site in SITES:
        for seed in SEEDS:
            tag, prefix, ok = run_site(site, seed)
            if not ok:
                print(f'FAIL site {site} seed {seed}', flush=True)
                continue
            med, frac, hard = diag(tag, prefix)
            runs.append(dict(site=site, seed=seed, max_p_median=float(med), frac_gt_half=float(frac), hard=hard.tolist()))
            print(f'ok site {site} seed {seed} maxp={med:.3f} f>.5={frac:.3f}', flush=True)
    # per-site seed agreement
    by_site = {}
    for r in runs:
        by_site.setdefault(r['site'], {})[r['seed']] = np.array(r['hard'])
    agree = []
    for site, d in sorted(by_site.items()):
        if len(d) < 2:
            continue
        a, b = d[SEEDS[0]], d[SEEDS[1]]
        ari = adjusted_rand_score(a, b)
        O = np.zeros((K, K))
        for x, y in zip(a, b):
            O[x, y] += 1
        ra, ca = linear_sum_assignment(-O)
        best = O[ra, ca].sum() / len(a)
        sa = np.bincount(a, minlength=K) / len(a)
        sb = np.bincount(b, minlength=K) / len(b)
        P = np.outer(sa, sb)
        rp, cp = linear_sum_assignment(-P)
        ch = P[rp, cp].sum()
        agree.append(dict(site=site, ari=ari, agree=float(best), chance=float(ch), excess=float(best - ch)))
        print(f'site {site}: ARI={ari:.3f} agree={best:.3f} chance={ch:.3f} excess={best-ch:+.3f}', flush=True)
    json.dump(dict(purpose='site-ladder tract-continuity probe', S=S, K=K, steps=STEPS, mg=MG,
                   marker_span=[lo, hi], sites=SITES, runs=runs, agreement=agree),
              open(OUTDIR / 'site_ladder.report.json', 'w'), indent=1, default=str)
    print('WROTE', OUTDIR / 'site_ladder.report.json')


if __name__ == '__main__':
    sys.exit(main())
