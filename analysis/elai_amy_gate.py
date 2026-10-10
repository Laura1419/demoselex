#!/usr/bin/env python3
"""Run the preregistered ELAI fits and fail-closed instrument gates; no biology."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import threading
import time

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/elai_amy'
BINARY = '/home/erikg/elai-fork/src/elai'
SEEDS = (20261010, 20261011)


def read_joint(path, n, k):
    p = np.full((n, k, k), np.nan)
    q = p.copy()
    counts = np.zeros(n, int)
    with path.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            i, a, b = (int(r[x]) for x in ('individual_index', 'k1', 'k2'))
            assert r['locus'] == 'virtual:1:104121000' and int(r['position']) == 104121000
            assert not np.isfinite(p[i, a, b]), 'duplicate joint row'
            p[i, a, b], q[i, a, b] = float(r['p_ordered']), float(r['q_unordered'])
            counts[i] += 1
    assert np.all(counts == k * k), 'incomplete individual joint'
    return p, q


def joint_errors(p, q):
    k = p.shape[1]
    return dict(ordered_sum_error=float(np.max(np.abs(p.sum((1, 2)) - 1))),
                q_symmetry_error=float(np.max(np.abs(q - q.transpose(0, 2, 1)))),
                ordered_symmetry_error=float(np.max(np.abs(p - p.transpose(0, 2, 1)))),
                q_triangle_sum_error=float(np.max(np.abs(q[:, np.triu_indices(k)[0], np.triu_indices(k)[1]].sum(1) - 1))),
                finite=bool(np.isfinite(p).all() and np.isfinite(q).all()),
                minimum=float(p.min()), maximum=float(p.max()))


def consume_full(fifo, n, k, m, result):
    """One-line-at-a-time FIFO drain: never materialize the full tensor."""
    counts = np.zeros(n, np.int64)
    with (OUT / 'samples.tsv').open() as f:
        samples = list(csv.DictReader(f, delimiter='\t'))
    with (OUT / 'markers.tsv').open() as f:
        markers = [r['id'] for r in csv.DictReader(f, delimiter='\t')]
    norm, error, rows = 0., 0., 0
    try:
        with fifo.open() as f:
            assert f.readline().startswith('individual_index\t')
            for line in f:
                fields = line.split('\t')
                i = int(fields[0]); a = int(fields[5]); b = int(fields[6])
                assert a == (rows % (k*k)) // k and b == rows % k
                if a == 0 and b == 0:
                    assert i == rows // (m*k*k)
                    # Upstream ignores BIMBAM IND names; numeric index is authoritative.
                    assert fields[1] in ('', samples[i]['sample'])
                    assert fields[2] == markers[(rows // (k*k)) % m]
                counts[i] += 1; rows += 1
                v = float(fields[7])
                assert np.isfinite(v) and 0 <= v <= 1
                norm += v
                if rows % (k*k) == 0:
                    error = max(error, abs(norm - 1)); norm = 0.
        result.update(rows=rows, individual_rows=counts.tolist(), norm_error=error,
                      complete=bool(rows == n*m*k*k and np.all(counts == m*k*k)))
    except Exception as exc:
        result.update(complete=False, error=str(exc), rows=rows)
        # The context above closes on a malformed stream; ELAI then fails loudly.
        return


def run_fit(work, k, seed, panel='panel', full=False, rem=None, steps=20, prefix=None):
    prefix = prefix or f'{panel}_k{k}_seed{seed}'
    log = OUT / f'{prefix}.stdout.txt'
    command = [BINARY, '-g', str(OUT / f'{panel}.inp'), '-p', '1', '-pos', str(OUT / f'{panel}.pos'),
               '-C', '1', '-c', str(k), '-w', '0', '-s', str(steps), '-mg', '10', '-R', str(seed),
               '-o', prefix]
    if steps:
        command += ['--pk1']
    command += ['--pkk-site', '1:104121000']
    if rem:
        command += ['-rem', str(rem)]
    n = json.loads((OUT / 'panel.json').read_text())['fitted_n']
    with (OUT / f'{panel}.inp').open() as f:
        f.readline(); m = int(f.readline())
    streamed = {}
    if full:
        command += ['--pkk']
        fifo = work / 'output' / f'{prefix}.pkk.tsv'
        os.mkfifo(fifo)
        thread = threading.Thread(target=consume_full, args=(fifo, n, k, m, streamed), daemon=True)
        thread.start()
    start = time.time()
    with log.open('w') as f:
        proc = subprocess.run(command, cwd=work, stdout=f, stderr=subprocess.STDOUT, timeout=900)
    if full:
        thread.join(timeout=30)
        assert not thread.is_alive(), 'FIFO did not finish'
    text = log.read_text()
    assert proc.returncode == 1 and '-elai: finished' in text, f'ELAI failed: {prefix}; stream={streamed}'
    source = work / 'output' / prefix
    p, q = read_joint(Path(str(source) + '.pkk-site.tsv'), n, k)
    # Assert retained panel, not only builder selection.
    retained = []
    for line in Path(str(source) + '.snpinfo.txt').read_text().splitlines():
        fields = line.split()
        if len(fields) >= 6:
            try:
                # Actual binary snpinfo header is rs,pos,chr,major,minor,maf,miss,n0,n1,n2.
                retained.append((fields[0], int(fields[2]), int(fields[1])))
            except ValueError:
                pass
    assert len(retained) == m
    excluded = sum(103915000 <= pos <= 104327000 for _, chrom, pos in retained)
    assert excluded == 0 and all(chrom == 1 for _, chrom, pos in retained)
    for suffix in ('.pkk-site.tsv', '.em.txt', '.log.txt', '.snpinfo.txt'):
        shutil.copyfile(Path(str(source) + suffix), OUT / (prefix + suffix))
    result = dict(prefix=prefix, k=k, seed=seed, panel=panel, command=command, cwd='DISPOSABLE_SCRATCH',
                  started_epoch=start, prereg_mtime_epoch=(ROOT / 'docs/elai-amy-prereg.md').stat().st_mtime,
                  elapsed_seconds=time.time()-start, returncode=proc.returncode, markers=m,
                  excluded_retained_markers=excluded, joint_errors=joint_errors(p, q), streamed=streamed)
    assert result['prereg_mtime_epoch'] < start, 'registration must precede fitting'
    if steps:
        comp = []
        with Path(str(source) + '.pk1.txt').open() as f:
            for line in f:
                x = np.fromstring(line, sep=' ').reshape(m, k) / 2
                assert np.isfinite(x).all()
                comp.append(x.mean(axis=0))
        comp = np.array(comp)
        assert comp.shape == (n, k)
        np.savetxt(OUT / (prefix + '.composition.tsv'), comp, delimiter='\t')
    (OUT / (prefix + '.run.json')).write_text(json.dumps(result, indent=2) + '\n')
    print(prefix, 'finished', round(result['elapsed_seconds'], 1), flush=True)
    return result


def silhouette(x, labels):
    return float(silhouette_score(x, labels)) if 1 < len(set(labels)) < len(labels) else None


def same_k(a, b, pa, pb):
    cost = np.abs(a[:, :, None] - b[:, None, :]).mean(0)
    rows, cols = linear_sum_assignment(cost)
    order = cols[np.argsort(rows)]
    b = b[:, order]
    pb = pb[:, order][:, :, order]
    return dict(ari=float(adjusted_rand_score(a.argmax(1), b.argmax(1))),
                mean_l1=float(np.abs(a-b).sum(1).mean()),
                virtual_argmax_agreement=float((pa.sum(2).argmax(1) == pb.sum(2).argmax(1)).mean()))


def gate_metrics(runs):
    panel = json.loads((OUT / 'panel.json').read_text())
    with (OUT / 'samples.tsv').open() as f:
        samples = list(csv.DictReader(f, delimiter='\t'))
    n = len(samples)
    pcs = np.loadtxt(OUT / 'independent_pcs.tsv', skiprows=1)
    primary = next(r for r in runs if r['k'] == 5 and r['seed'] == SEEDS[0] and r['panel'] == 'panel')
    comp, joints = {}, {}
    for r in runs:
        if r['prefix'].startswith('frozen'):
            continue
        key = (r['k'], r['seed'], r['panel'])
        comp[key] = np.loadtxt(OUT / (r['prefix'] + '.composition.tsv'))
        joints[key] = read_joint(OUT / (r['prefix'] + '.pkk-site.tsv'), n, r['k'])[0]
    a = comp[5, SEEDS[0], 'panel']; labels = a.argmax(1)
    groups = [s['group'] for s in samples]; countries = [s['country'] for s in samples]
    date = np.array([float(s['BP']) for s in samples])
    date_fit = np.column_stack((np.ones(n), a)) @ np.linalg.lstsq(np.column_stack((np.ones(n), a)), date, rcond=None)[0]
    known = dict(nonempty_assignments=len(set(labels)), silhouette=silhouette(a, labels),
                 pc_silhouette=silhouette(pcs, labels), group_nmi=float(normalized_mutual_info_score(groups, labels)),
                 country_nmi=float(normalized_mutual_info_score(countries, labels)),
                 date_bin_nmi=float(normalized_mutual_info_score(np.floor(date/2000), labels)),
                 date_composition_r2=float(1-((date-date_fit)**2).sum()/((date-date.mean())**2).sum()),
                 group_pc_silhouette=silhouette(pcs, groups), country_pc_silhouette=silhouette(pcs, countries),
                 group_composition_silhouette=silhouette(a, groups), country_composition_silhouette=silhouette(a, countries))
    known['pass'] = all(known[x] is not None and known[x] > threshold for x, threshold in
                        [('silhouette', .10), ('pc_silhouette', .10), ('group_nmi', .10), ('country_nmi', .05)])
    stability = dict(seed=[], cross_k=[], thinning=[])
    for k in (3, 5, 8):
        x, y = (k, SEEDS[0], 'panel'), (k, SEEDS[1], 'panel')
        if x in comp and y in comp:
            v = same_k(comp[x], comp[y], joints[x], joints[y]); v['k'] = k
            v['pass'] = v['ari'] >= .8 and v['mean_l1'] <= .2 and v['virtual_argmax_agreement'] >= .8
            stability['seed'].append(v)
    for seed in SEEDS:
        for k, l in ((3, 5), (5, 8)):
            x, y = (k, seed, 'panel'), (l, seed, 'panel')
            if x in comp and y in comp:
                c, d = comp[x], comp[y]
                # Normalize off-diagonal co-composition similarities to remove K-dependent scale.
                g, h = c @ c.T, d @ d.T
                g /= np.linalg.norm(g); h /= np.linalg.norm(h)
                v = dict(seed=seed, k_pair=[k, l], ari=float(adjusted_rand_score(c.argmax(1), d.argmax(1))),
                         gram_difference=float(np.linalg.norm(g-h)))
                v['pass'] = v['ari'] >= .6 and v['gram_difference'] <= .3
                stability['cross_k'].append(v)
    x, y = (5, SEEDS[0], 'panel'), (5, SEEDS[0], 'thin')
    if y in comp:
        v = same_k(comp[x], comp[y], joints[x], joints[y])
        v['pass'] = v['ari'] >= .8 and v['mean_l1'] <= .2 and v['virtual_argmax_agreement'] >= .8
        stability['thinning'].append(v)
    stability['comparisons_completed'] = sum(len(v) for v in stability.values() if isinstance(v, list))
    stability['comparisons_required'] = 8
    stability['pass'] = stability['comparisons_completed'] == 8 and all(v['pass'] for name in ('seed', 'cross_k', 'thinning') for v in stability[name])
    coverage = dict(fitted_n=n, usable_individuals=panel['usable_individuals'],
                    usable_fraction=panel['usable_individuals']/n, dated_virtual_output=n,
                    full_export_rows=primary['streamed'].get('rows', 0),
                    full_export_complete=primary['streamed'].get('complete', False))
    coverage['pass'] = coverage['usable_fraction'] >= .95 and n >= 100 and coverage['full_export_complete']
    errors = [r['joint_errors'] for r in runs]
    sanity = {name: max(e[name] for e in errors) for name in ('ordered_sum_error', 'q_symmetry_error', 'ordered_symmetry_error', 'q_triangle_sum_error')}
    linked = OUT / 'frozen_linked.pkk-site.tsv'; unlinked = OUT / 'frozen_unlinked.pkk-site.tsv'
    completed_prefixes = {r['prefix'] for r in runs}
    sanity['bridge_control_completed'] = int({'frozen_linked', 'frozen_unlinked'} <= completed_prefixes)
    sanity['bridge_max_change'] = None
    if sanity['bridge_control_completed']:
        p, _ = read_joint(linked, n, 5); q, _ = read_joint(unlinked, n, 5)
        sanity['bridge_max_change'] = float(np.abs(p-q).max())
    sanity['pass'] = (sanity['ordered_sum_error'] <= 1e-10 and sanity['q_symmetry_error'] <= 1e-12
                      and sanity['q_triangle_sum_error'] <= 1e-10 and all(e['finite'] and e['minimum'] >= 0 and e['maximum'] <= 1 for e in errors)
                      and sanity['bridge_max_change'] is not None and sanity['bridge_max_change'] > 1e-6)
    exclusion = dict(builder_core_markers=panel['excluded_markers'], retained_core_markers=sum(r['excluded_retained_markers'] for r in runs),
                     site_inside_core=panel['site_inside_core'], left_flank=panel['left_flank'], right_flank=panel['right_flank'],
                     pc_overlap=panel['pc_overlap'])
    exclusion['pass'] = exclusion['builder_core_markers'] == exclusion['retained_core_markers'] == exclusion['pc_overlap'] == 0 and exclusion['site_inside_core'] and exclusion['left_flank'] < 104121000 < exclusion['right_flank']
    gates = dict(known_structure=known, stability=stability, coverage=coverage, virtual_sanity=sanity, exclusion=exclusion)
    result = dict(gates=gates, verdict='PASS' if all(g['pass'] for g in gates.values()) else 'STOP',
                  runs=len(runs), execution_complete=len(runs) == 9, disk_free_bytes=shutil.disk_usage('/').free,
                  prereg_sha256=hashlib.sha256((ROOT / 'docs/elai-amy-prereg.md').read_bytes()).hexdigest())
    (OUT / 'gate_metrics.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(result, indent=2), flush=True)
    return result


def unlinked_parameters(source, destination, right):
    """Change only gap rk, preserving the upstream line-oriented parameter format."""
    lines = source.read_text().splitlines(keepends=True)
    rk = lines[2].split()
    original = rk[right]
    rk[right] = '1.0'
    lines[2] = ' '.join(rk) + '\n'
    destination.write_text(''.join(lines))
    return original


def run_controls(work):
    primary = OUT / 'panel_k5_seed20261010.em.txt'
    rows = np.loadtxt(OUT / 'markers.tsv', delimiter='\t', skiprows=1, usecols=[2], dtype=int)
    right = int(np.searchsorted(rows, 104121000))
    unlinked = OUT / 'unlinked_control.em.txt'
    original = unlinked_parameters(primary, unlinked, right)
    runs = [run_fit(work, 5, SEEDS[0], rem=primary, steps=0, prefix='frozen_linked'),
            run_fit(work, 5, SEEDS[0], rem=unlinked, steps=0, prefix='frozen_unlinked')]
    (OUT / 'bridge_control.json').write_text(json.dumps(dict(rk_token_index=2*len(rows)+right, right_marker_index=right,
                                                    original_rk=original, unlinked_rk=1.0), indent=2)+'\n')
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--metrics-only', action='store_true')
    ap.add_argument('--controls-only', action='store_true')
    args = ap.parse_args()
    if args.metrics_only:
        gate_metrics([json.loads(p.read_text()) for p in OUT.glob('*.run.json')]); return
    started = time.time()
    # Invalidate stale metrics before any fit; an aborted primary can never reuse a prior pass.
    pending = dict(verdict='STOP', runs=0, execution_complete=False,
                   gates={name: {'pass': False, 'tests_completed': 0, 'reason': 'execution incomplete'}
                          for name in ('known_structure', 'stability', 'coverage', 'virtual_sanity', 'exclusion')})
    (OUT / 'gate_metrics.json').write_text(json.dumps(pending, indent=2) + '\n')
    work = Path(tempfile.mkdtemp(prefix='elai-amy-', dir='/tmp'))
    (work / 'output').mkdir()
    runs = []
    try:
        if args.controls_only:
            runs = [json.loads(p.read_text()) for p in OUT.glob('*.run.json') if not p.name.startswith('frozen')]
        else:
            runs.append(run_fit(work, 5, SEEDS[0], full=True))
            remaining = [(k, seed, 'panel') for k in (3, 5, 8) for seed in SEEDS if (k, seed) != (5, SEEDS[0])]
            remaining.append((5, SEEDS[0], 'thin'))
            with ThreadPoolExecutor(max_workers=3) as pool:
                futures = [pool.submit(run_fit, work, k, seed, panel) for k, seed, panel in remaining]
                for future in futures:
                    runs.append(future.result())
        runs.extend(run_controls(work))
    except Exception as exc:
        (OUT / 'execution_failure.json').write_text(json.dumps(dict(error=str(exc), elapsed_seconds=time.time()-started), indent=2)+'\n')
        print('FAIL-CLOSED:', repr(exc), flush=True)
    finally:
        shutil.rmtree(work)
        (OUT / 'scratch_cleanup.json').write_text(json.dumps(dict(path=str(work), removed=not work.exists(), elapsed_seconds=time.time()-started), indent=2)+'\n')
    if runs:
        gate_metrics(runs)


if __name__ == '__main__':
    main()
