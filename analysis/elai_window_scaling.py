"""Diagnostic probe: does the local-state posterior depend on window scale?

Hypothesis under test (raised 2026-10-10): the dense panel spans a 35 Mb region
(chr1:85-120 Mb) while ancestry tract lengths at Neolithic-age admixture are
~0.5-1 Mb. If the switch/recombination model is effectively dead at that scale
(-mg measured inert at 5/10/100/1000), the posterior at the AMY virtual site may
be flat because the model is assigning one block-wide state across 35 Mb, not
because the states are fundamentally unidentifiable.

Design: subset the SAME retained dense markers and SAME 1,000 dated individuals
to windows centred on the AMY virtual site (chr1:104,121,000), run ELAI with the
settled configuration (K=4, S=1, 100 steps) and the --pkk-site export only (no
full marginal export), several seeds per width, and measure posterior
concentration and label-permutation-invariant seed agreement.

This is a diagnostic probe, not a pre-registered gate. It reuses the dense
build's inputs read-only and writes only under data/derived/elai_ws/.
"""
import argparse, json, subprocess, sys
from bisect import bisect_left, bisect_right
from itertools import permutations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

BINARY = '/home/erikg/elai-fork/src/elai'
SRC = Path('data/derived/elai_amy_dense/work_long')
CENTER = 104121000
SITE = '1:104121000'
OUT = Path('data/derived/elai_ws')


def build_window(half_bp, tag, reuse=True):
    dst = OUT / tag
    dst.mkdir(parents=True, exist_ok=True)
    if reuse and all((dst / f).exists() for f in ('panel.inp', 'panel.pos', 'genotypes.txt')):
        return None
    markers = json.load(open(SRC / 'markers.json'))
    samples = json.load(open(SRC / 'samples.json'))
    positions = [m[2] for m in markers]
    lo, hi = CENTER - half_bp, CENTER + half_bp
    a, b = bisect_left(positions, lo), bisect_right(positions, hi)
    subset = markers[a:b]
    n_left = sum(1 for m in subset if m[2] < CENTER - 206000)
    n_right = sum(1 for m in subset if m[2] > CENTER + 206000)
    if len(subset) < 200 or n_left < 50 or n_right < 50:
        raise RuntimeError(f'{tag}: insufficient flanking coverage left={n_left} right={n_right} total={len(subset)}')
    with open(SRC / 'genotypes.txt') as g:
        lines = g.readlines()
    assert len(lines) == len(markers), 'genotype rows must match marker list'
    keep = set(range(a, b))
    with (dst / 'genotypes.txt').open('w') as t:
        for i in sorted(keep):
            t.write(lines[i])
    with (dst / 'panel.inp').open('w') as f:
        f.write(f"{len(samples)}\n{len(subset)}\nIND " + ' '.join(s['sample'] for s in samples) + '\n')
        for i in sorted(keep):
            f.write(lines[i])
    with (dst / 'panel.pos').open('w') as f:
        for m in subset:
            f.write(f'{m[3]} {m[2]} {m[1]}\n')
    return dict(tag=tag, half_bp=half_bp, width_mb=2 * half_bp / 1e6, markers=len(subset),
                flank_left=n_left, flank_right=n_right, n=len(samples))


def run_fit(tag, seed, k=4, steps=100, mg=10, S=1):
    dst = OUT / tag
    prefix = f'w_s{S}_k{k}_seed{seed}_mg{mg}'
    cmd = [BINARY, '-g', 'panel.inp', '-p', '1', '-pos', 'panel.pos', '-C', str(S), '-c', str(k),
           '-w', '0', '-s', str(steps), '-mg', str(mg), '-R', str(seed), '-o', prefix,
           '--pkk-site', SITE]
    with (dst / f'{prefix}.launch.log').open('w') as log:
        proc = subprocess.run(cmd, cwd=dst, stdout=log, stderr=subprocess.STDOUT, timeout=7200)
    log_text = (dst / f'{prefix}.launch.log').read_text()
    ok = proc.returncode == 1 and '-elai: finished' in log_text
    trace = []
    for ln in log_text.splitlines():
        if ln.startswith(' ### Exact update'):
            parts = ln.split()
            trace.append((int(parts[4]), float(parts[5])))
    final = trace[-1][1] if trace else None
    return dict(prefix=prefix, returncode=proc.returncode, finished=ok, cmd=' '.join(cmd),
                final_likelihood=final, steps_logged=len(trace))


def diagnostics(tag, prefix, k=4):
    f = OUT / tag / 'output' / f'{prefix}.pkk-site.tsv'
    if not f.exists():
        raise RuntimeError(f'missing {f}')
    d = pd.read_csv(f, sep='\t')
    assert abs(d.groupby('individual_index').p_ordered.sum().median() - 1.0) < 1e-6
    per = d.groupby('individual_index').p_ordered
    mx = per.max()
    ent = per.apply(lambda p: -(p * np.log(p)).sum())
    marg = d.pivot_table(index='individual_index', columns='k1', values='p_ordered', aggfunc='sum')
    hard = marg.idxmax(axis=1).to_numpy()
    return dict(max_p_median=float(mx.median()), frac_gt_half=float((mx > .5).mean()),
                frac_gt_third=float((mx > .3).mean()), entropy_median=float(ent.median()),
                uniform_entropy=float(np.log(k * k)), hard=hard)


def agreement(diags, k):
    """Label-permutation-invariant pairwise agreement via Hungarian alignment
    (brute force is infeasible for K>8)."""
    from scipy.optimize import linear_sum_assignment
    pairs = []
    for i in range(len(diags)):
        for j in range(i + 1, len(diags)):
            a, b = diags[i]['hard'], diags[j]['hard']
            ari = adjusted_rand_score(a, b)
            O = np.zeros((k, k))
            for x, y in zip(a, b):
                O[x, y] += 1
            ra, ca = linear_sum_assignment(-O)
            best = float(O[ra, ca].sum() / len(a))
            sa = np.bincount(a, minlength=k) / len(a)
            sb = np.bincount(b, minlength=k) / len(b)
            P = np.outer(sa, sb)
            rp, cp = linear_sum_assignment(-P)
            chance = float(P[rp, cp].sum())
            pairs.append(dict(ari=ari, best_permuted_agreement=best, chance_agreement=chance,
                              excess=best - chance, perm=list(int(c) for c in ca)))
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--widths', default='0.5,1,2,5', help='full window widths in Mb')
    ap.add_argument('--seeds', default='20261020,20261021,20261022')
    ap.add_argument('--k', type=int, default=4)
    ap.add_argument('--steps', type=int, default=100)
    ap.add_argument('--mg-sweep', default='', help='comma list of extra -mg values to probe at the first width')
    ap.add_argument('--report-suffix', default='', help='suffix for the report json (allows parallel launches)')
    ap.add_argument('--S', type=int, default=1, help='upper-layer source clusters (-C); 1 = current default')
    args = ap.parse_args()
    suffix = args.report_suffix
    seeds = [int(s) for s in args.seeds.split(',')]
    report = dict(purpose='window-scale diagnostic probe (not a prereg gate)', k=args.k, S=args.S, steps=args.steps,
                  center=CENTER, windows=[], mg_sensitivity=[])
    first_tag = None
    for w in [float(x) for x in args.widths.split(',')]:
        tag = f'w{int(w * 1000)}kb'
        if first_tag is None:
            first_tag = tag
        meta = build_window(int(w * 1e6 / 2), tag)
        runs, diags = [], []
        for seed in seeds:
            r = run_fit(tag, seed, k=args.k, steps=args.steps, S=args.S)
            runs.append(r)
            if r['finished']:
                diags.append(diagnostics(tag, r['prefix'], k=args.k))
        rep = dict(**(meta or dict(tag=tag, half_bp=int(w * 1e6 / 2), width_mb=w, markers=None,
                                   flank_left=None, flank_right=None, n=1000, reused=True)),
                   runs=runs, seed_agreement=agreement(diags, args.k) if len(diags) > 1 else [])
        for d in diags:
            d.pop('hard')
        rep['diagnostics'] = diags
        report['windows'].append(rep)
        print(json.dumps(rep), flush=True)
    for mgv in [int(x) for x in args.mg_sweep.split(',') if x]:
        tag = first_tag
        r = run_fit(tag, seeds[0], k=args.k, steps=args.steps, mg=mgv, S=args.S)
        d = diagnostics(tag, r['prefix'], k=args.k) if r['finished'] else None
        d.pop('hard', None)
        report['mg_sensitivity'].append(dict(tag=tag, mg=mgv, run=r, diagnostics=d))
        print(json.dumps(report['mg_sensitivity'][-1]), flush=True)
    json.dump(report, open(OUT / f'window_scaling.report{suffix}.json', 'w'), indent=1, default=str)
    print('WROTE', OUT / f'window_scaling.report{suffix}.json')


if __name__ == '__main__':
    sys.exit(main())
