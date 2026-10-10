#!/usr/bin/env python3
"""Persistent dense ELAI driver; validation fixtures never authorize biology."""
import argparse
import contextlib
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

# Fix numerical threading before numpy/sklearn imports, including direct CLI use.
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
sys.dont_write_bytecode = True
import numpy as np
import elai_amy_dense as dense

BASE = dense.OUT
DEFAULT_WORK = BASE / 'work_long'


def save(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


@contextlib.contextmanager
def lock(path, exclusive, blocking=True):
    with path.open('a') as f:
        flags = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
        if not blocking:
            flags |= fcntl.LOCK_NB
        try:
            fcntl.flock(f, flags)
        except BlockingIOError:
            raise RuntimeError(f'already running: {path}') from None
        yield


def load(work):
    panel = json.loads((work / 'panel.json').read_text())
    if hashlib.sha256(dense.PREREG.read_bytes()).hexdigest() != panel['prereg_sha256']:
        raise RuntimeError('preregistration differs from built input')
    samples = json.loads((work / 'samples.json').read_text())
    markers = json.loads((work / 'markers.json').read_text())
    pc = np.load(work / 'independent_pcs.npy', allow_pickle=False)
    assert len(samples) == panel['fitted_n'] and len(markers) == panel['retained_markers']
    assert pc.shape == (len(samples), 5)
    return panel, samples, markers, pc


def fixture_build(work):
    """Deterministic real-input subset: 12 individuals, 80 spaced markers + flanks."""
    panel, samples, markers, pc = load(DEFAULT_WORK)
    indices = set(np.linspace(0, len(markers)-1, 80, dtype=int).tolist())
    right = int(np.searchsorted([m[2] for m in markers], 104121000))
    indices.update((right-1, right))
    selected = sorted(indices)
    subset = [markers[i] for i in selected]
    samples = samples[:12]
    with (DEFAULT_WORK / 'genotypes.txt').open() as source, (work / 'genotypes.txt').open('w') as target:
        wanted = set(selected)
        for i, line in enumerate(source):
            if i in wanted:
                target.write(' '.join(line.split()[:len(samples)+1]) + '\n')
    with (work / 'panel.inp').open('w') as f, (work / 'genotypes.txt').open() as g:
        f.write(f'{len(samples)}\n{len(subset)}\nIND ' + ' '.join(s['sample'] for s in samples) + '\n')
        shutil.copyfileobj(g, f)
    with (work / 'panel.pos').open('w') as f:
        for m in subset:
            f.write(f'{m[3]} {m[2]} {m[1]}\n')
    panel.update(fitted_n=len(samples), dated_n=len(samples), retained_markers=len(subset),
                 individual_call_rates=panel['individual_call_rates'][:12], usable_individuals=12,
                 validation_fixture=True, fixture_source=str(DEFAULT_WORK),
                 input_bytes=(work / 'panel.inp').stat().st_size)
    dense.save('panel.json', panel)
    return panel, samples, subset, pc[:12]


def build(work, fixture):
    if (work / 'build.summary.json').exists() or (work / 'panel.inp').exists():
        raise RuntimeError('build already exists; use it read-only or choose a new --work')
    (work / 'output').mkdir(exist_ok=True)
    panel, samples, markers, pc = fixture_build(work) if fixture else dense.build(work)
    # The original builder deletes its intermediate. Reconstitute exactly its
    # genotype lines, WITHOUT the three BIMBAM header lines, and retain them.
    if not fixture:
        with (work / 'panel.inp').open() as f, (work / 'genotypes.txt').open('w') as g:
            for _ in range(3):
                next(f)
            shutil.copyfileobj(f, g, length=1024*1024)
        assert panel['retained_markers'] >= 50000 and panel['pc_markers'] >= 5000 and panel['fitted_n'] >= 1000
    save(work / 'samples.json', samples)
    save(work / 'markers.json', markers)
    np.save(work / 'independent_pcs.npy', pc, allow_pickle=False)
    return dict(retained_markers=len(markers), fitted_n=len(samples), pc_markers=panel['pc_markers'],
                validation_fixture=fixture, sizes_bytes={name: (work / name).stat().st_size
                    for name in ('panel.inp', 'panel.pos', 'genotypes.txt')})


def run_fit(work, args):
    panel, samples, markers, pc = load(work)
    prefix = f'panel_k{args.k}_seed{args.seed}'
    with lock(work / (prefix + '.lock'), True, False):
        if list(work.glob(prefix + '.*.json')) or list((work / 'output').glob(prefix + '*')):
            raise RuntimeError(f'prefix already used: {prefix}; preserve it and choose another work/seed')
        dense.STEPS = args.steps  # Process-local override: never edit the prereg or source constant.
        dense.seen_steps.clear()
        result = dense.fit(work, panel, markers, args.k, args.seed, args.wall,
                           full_export=args.export == 'joint')
        result['export_mode'] = args.export
        result['validation_fixture'] = bool(panel.get('validation_fixture'))
        result['requested_steps'] = args.steps
        result['cwd'] = str(work)
        result['streamed_markers'] = len(markers)
        result['streamed_individuals'] = panel['fitted_n'] if (
            result['full_export_complete'] or result.get('marginal_export', {}).get('complete')) else None
        save(work / (prefix + '.run.json'), result)
        if not result['complete'] or not (result['full_export_complete'] or result.get('marginal_export', {}).get('complete')):
            raise RuntimeError(f'fit/export failed: inspect {prefix}.run.json')
    return result


def status(passed):
    return dict(status='PASS' if passed else 'FAIL', pass_gate=bool(passed))


def aggregate(parts):
    if any(p['status'] == 'FAIL' for p in parts):
        return 'FAIL'
    return 'PASS' if all(p['status'] == 'PASS' for p in parts) else 'PARTIAL'


def gates(work):
    panel, samples, markers, pc = load(work)
    keys = ((5, 20261010), (5, 20261011), (3, 20261010))
    runs = []
    for k, seed in keys:
        path = work / f'panel_k{k}_seed{seed}.run.json'
        if path.exists():
            runs.append(json.loads(path.read_text()))
        else:
            runs.append(dict(prefix=f'panel_k{k}_seed{seed}', k=k, seed=seed,
                             complete=False, full_export_complete=False,
                             full_export_rows=0,
                             convergence=dict(status='NOT TESTED', pass_gate=False,
                                              completed_steps=0, trace=[])))
    # Reuse the preregistered metric implementation. Its validity flag meant
    # complete joint export. Marginals supply compositions, but this adapter
    # NEVER claims joint coverage and NEVER accepts nonconverged fit labels.
    adapted = []
    for r in runs:
        a = dict(r)
        a['full_export_complete'] = bool(r.get('convergence', {}).get('pass_gate') and
            (r.get('full_export_complete') or r.get('marginal_export', {}).get('complete')))
        adapted.append(a)
    result = dense.metrics(panel, samples, pc, adapted)
    g = result['gates']
    primary = runs[0]
    missing = dict(status='NOT TESTED', pass_gate=False)
    converged = primary.get('convergence', {}).get('pass_gate', False)
    known = g['known_structure']
    if known['status'] != 'NOT TESTED':
        # Exact prereg inequalities; dense.metrics originally used >= for NMI.
        known.update(status(known['nonempty_assignments'] >= 2 and
            known['pc_silhouette'] is not None and known['pc_silhouette'] >= .1 and
            known['composition_silhouette'] is not None and known['composition_silhouette'] >= .1 and
            known['group_nmi'] > .1 and known['country_nmi'] > .05))
    known['composition_source'] = primary.get('export_mode')
    known['singleton_robust_measure'] = 'country adjusted mutual information (diagnostic; no prereg threshold)'
    known['reason_if_untested'] = None if converged else 'No converged primary: labels not accepted'
    stab = g['stability']
    stab['status'] = aggregate([stab[x] for x in ('seed', 'cross_k', 'virtual_argmax', 'thinning_refit')])
    stab['pass_gate'] = stab['status'] == 'PASS'
    coverage = g['coverage']
    stream_complete = primary.get('full_export_complete', False)
    marginal_complete = primary.get('marginal_export', {}).get('complete', False)
    coverage.update(full_export_rows=primary.get('full_export_rows', 0),
                    dated_complete_export=panel['dated_n'] if stream_complete else 0)
    coverage['components'] = dict(
        usable_input=dict(status(panel['usable_individuals']/panel['fitted_n'] >= .95),
                          usable_individuals=panel['usable_individuals']),
        dated_complete_joint=dict(status(panel['dated_n'] >= 100), dated_n=panel['dated_n']) if stream_complete else
            dict(missing, reason='complete joint export unavailable'),
        full_ordered_export=dict(status(stream_complete), rows=primary.get('full_export_rows', 0)) if primary.get('export_mode') == 'joint' else
            dict(missing, reason='Export-size-driven omission: marginal mode is not N*M*K^2 joint coverage'),
        marginal_export=dict(status(marginal_complete), **primary.get('marginal_export', {})) if primary.get('export_mode') == 'marginal' else dict(missing))
    required = [coverage['components'][x] for x in ('usable_input', 'dated_complete_joint', 'full_ordered_export')]
    coverage['status'] = aggregate(required)
    coverage['pass_gate'] = coverage['status'] == 'PASS'
    # Arithmetic is observable even for a finished but nonconverged smoke fit.
    virtual = g['virtual_site']
    errors = [r['joint_errors'] for r in runs if r.get('joint_errors')]
    arithmetic = dict(missing)
    if errors:
        maxima = {x: max(e[x] for e in errors) for x in
                  ('ordered_sum_error', 'q_symmetry_error', 'q_triangle_sum_error')}
        arithmetic = dict(status(maxima['ordered_sum_error'] <= 1e-10 and maxima['q_symmetry_error'] <= 1e-12 and
            maxima['q_triangle_sum_error'] <= 1e-10 and all(e['finite'] and e['minimum'] >= 0 and e['maximum'] <= 1 for e in errors)), **maxima)
        virtual.update(maxima)
    virtual['components'] = dict(arithmetic=arithmetic, integration=dict(missing, reason='Frozen linkage replays not run by this driver'))
    virtual['status'] = aggregate(list(virtual['components'].values()))
    virtual['pass_gate'] = virtual['status'] == 'PASS'
    # Retained SNP-info is written before EM; audit every ID/coordinate, not just counts.
    checked = 0
    for r in runs:
        p = work / 'output' / (r['prefix'] + '.snpinfo.txt')
        if p.exists():
            retained = [(f[0], int(f[1]), int(f[2])) for line in p.read_text().splitlines()[1:] if (f := line.split())]
            assert retained == [(m[3], m[2], m[1]) for m in markers], 'retained marker mismatch'
            checked += len(retained)
    exclusion = g['exclusion']
    exclusion['retained_markers_checked'] = checked
    if checked:
        exclusion['retained_core_markers'] = 0
        exclusion.update(status(panel['excluded_markers'] == 0 and panel['pc_overlap'] == 0 and
            all(m[1] == 1 and not dense.excluded(m[2]) for m in markers) and
            panel['left_flank'] < 104121000 < panel['right_flank'] and dense.excluded(104121000)))
    g['em_convergence'] = primary.get('convergence', dict(missing))
    result['verdict'] = aggregate(list(g.values()))
    result['validation_fixture'] = bool(panel.get('validation_fixture'))
    result['biological_use_authorized'] = result['verdict'] == 'PASS' and not panel.get('validation_fixture', False)
    result['fits'] = sum(r.get('complete', False) for r in runs)
    save(work / 'gate_metrics.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('build', 'fit', 'gates'), required=True)
    parser.add_argument('--work', type=Path, default=DEFAULT_WORK)
    parser.add_argument('--fixture', action='store_true', help='build a deterministic small subset of work_long')
    parser.add_argument('--k', type=int, choices=(3, 5))
    parser.add_argument('--seed', type=int)
    parser.add_argument('--steps', type=int, default=200)
    parser.add_argument('--wall', type=float, default=15000)
    parser.add_argument('--export', choices=('marginal', 'joint'), default='marginal')
    args = parser.parse_args()
    if args.mode == 'fit' and (args.k is None or args.seed is None or args.steps < 1 or args.wall <= 0):
        parser.error('fit needs --k, --seed, positive --steps and --wall')
    if args.fixture and args.mode != 'build':
        parser.error('--fixture is only for build')
    work = args.work.resolve()
    if not work.is_relative_to(BASE.resolve()) or work == BASE.resolve():
        parser.error('--work must be a subdirectory of data/derived/elai_amy_dense/')
    if args.fixture and work == DEFAULT_WORK:
        parser.error('fixture must use its own --work, never work_long')
    if args.mode != 'build' and not (work / 'build.summary.json').exists():
        parser.error('complete --mode build first')
    work.mkdir(parents=True, exist_ok=True)
    dense.OUT = work
    sys.path.insert(0, str(dense.ROOT))
    import scripts.provenance as provenance
    provenance.get_git_info = lambda: {'disabled': 'operator prohibits all git commands'}
    prefix = f'panel_k{args.k}_seed{args.seed}' if args.mode == 'fit' else args.mode
    summary_path = work / (prefix + '.summary.json')
    if summary_path.exists() and args.mode != 'gates':
        parser.error(f'preserving previous step summary: {summary_path}; choose a new work/seed')
    recorder = provenance.RunRecorder(f'elai-amy-long-{work.name}-{prefix}-{time.time_ns()}')
    recorder.set_seed(args.seed if args.seed is not None else dense.SEED)
    recorder.record_parameter('arguments', {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()})
    for p in (Path(__file__), Path(dense.__file__), dense.PREREG, Path(dense.BINARY)):
        recorder.record_input(p)
    if args.mode == 'build' and not args.fixture:
        recorder.record_parameter('raw_input_stats', {str(p): dict(size=p.stat().st_size, mtime=p.stat().st_mtime) for p in dense.paths()})
    else:
        for name in ('panel.json', 'samples.json', 'markers.json', 'independent_pcs.npy'):
            if (work / name).exists():
                recorder.record_input(work / name)
    start = time.time()
    summary_path = work / (prefix + '.summary.json')
    summary = dict(mode=args.mode, work=str(work), started_epoch=start)
    try:
        with lock(work / '.build.lock', args.mode == 'build'):
            with lock(work / '.gates.lock', True) if args.mode == 'gates' else contextlib.nullcontext():
                fn = {'build': lambda: build(work, args.fixture), 'fit': lambda: run_fit(work, args), 'gates': lambda: gates(work)}[args.mode]
                summary['result'] = fn()
        summary['success'] = True
    except Exception as exc:
        summary.update(success=False, error=repr(exc))
        raise
    finally:
        summary['elapsed_seconds'] = time.time()-start
        save(summary_path, summary)
        recorder.record_output(summary_path)
        recorder.record_parameter('input_sizes_bytes', {name: (work / name).stat().st_size for name in ('panel.inp', 'panel.pos', 'genotypes.txt') if (work / name).exists()})
        recorder.save()
        print('[observation]', json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
