#!/usr/bin/env python3
"""Trace-only ELAI settling diagnostic; never authorizes biological use."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'data/derived/elai_amy_dense'
BINARY = Path('/home/erikg/elai-fork/src/elai')
# name, S, K, warm steps, mixture generations, exact steps
CONFIGS = (
    ('baseline', 1, 5, 0, 10, 70),
    ('s1_k3', 1, 3, 0, 10, 100),
    ('s1_k4', 1, 4, 0, 10, 90),
    ('s2_k2', 2, 2, 0, 10, 90),
    ('s2_k3', 2, 3, 0, 10, 50),
    ('s3_k2', 3, 2, 0, 10, 50),
    ('mg5', 1, 5, 0, 5, 70),
    ('mg100', 1, 5, 0, 100, 70),
    ('mg1000', 1, 5, 0, 1000, 70),
    ('warm1', 1, 5, 1, 10, 70),
)
SEED = 20261010
TOLERANCE = 1e-4
DECAY_RATIO = 0.5
NEGATIVE_LIMIT = 0.05
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
EXACT = re.compile(r'^\s*### Exact update\s+(\d+)\s+(\d+)\s+(' + NUMBER + r')\s+(' + NUMBER + r')\s*$')


def parse_trace(text: str) -> list[dict]:
    trace = []
    for line in text.splitlines():
        match = EXACT.match(line)
        if not match:
            if '### Exact update' in line:
                raise ValueError(f'malformed exact update: {line}')
            continue
        run, index, likelihood, increment = match.groups()
        if int(run) != 0 or int(index) != len(trace):
            raise ValueError('expected one EM run with contiguous zero-based exact steps')
        likelihood, increment = float(likelihood), float(increment)
        if not math.isfinite(likelihood) or not math.isfinite(increment):
            raise ValueError('nonfinite trace')
        previous = trace[-1]['likelihood'] if trace else None
        # The first reported increment is L_0 - 0, not an EM improvement.
        delta = likelihood - previous if previous is not None else None
        trace.append(dict(step=int(index) + 1, elai_step=int(index), likelihood=likelihood,
                          reported_increment=increment, absolute_increment=delta,
                          relative_increment=delta / max(1.0, abs(previous)) if previous is not None else None))
    return trace


def decision(trace: list[dict], complete: bool) -> dict:
    increments = trace[1:]
    quarter = len(increments) // 4
    second = increments[len(increments) // 2:]
    first_mean = statistics.mean(abs(t['relative_increment']) for t in increments[:quarter]) if quarter else None
    last_mean = statistics.mean(abs(t['relative_increment']) for t in increments[-quarter:]) if quarter else None
    ratio = last_mean / first_mean if first_mean else None
    decays = bool(first_mean is not None and first_mean > 0 and last_mean <= DECAY_RATIO * first_mean)
    negative = statistics.mean(t['relative_increment'] < 0 for t in second) if second else None
    hits = [trace[i]['step'] for i in range(3, len(trace))
            if all(0 <= t['relative_increment'] < TOLERANCE for t in trace[i-2:i+1])]
    final_window = bool(hits and hits[-1] == trace[-1]['step'])
    return dict(completed_steps=len(trace), usable_increments=len(increments), quarter_size=quarter,
                first_quarter_mean_abs_relative_increment=first_mean,
                last_quarter_mean_abs_relative_increment=last_mean,
                last_over_first_ratio=ratio, decays=decays,
                second_half_mean_abs_increment=statistics.mean(abs(t['absolute_increment']) for t in second) if second else None,
                second_half_mean_abs_relative_increment=statistics.mean(abs(t['relative_increment']) for t in second) if second else None,
                second_half_negative_fraction=negative,
                first_qualifying_step=hits[0] if hits else None,
                final_three_window_met=final_window, complete=complete,
                settles=bool(complete and decays and negative is not None and negative < NEGATIVE_LIMIT and final_window))


def render(result: dict) -> str:
    rows = ['[observation] ELAI settling sweep trace report (not biological authorization).']
    for config in result['configurations']:
        stats = config['statistics']
        rows.extend(['', f"{config['name']}: completed={stats['completed_steps']}/{config['requested_steps']} "
                     f"complete={stats['complete']} elapsed_s={config['elapsed_seconds']} "
                     f"returncode={config['returncode']}",
                     '[observation] ' + json.dumps(stats, sort_keys=True, allow_nan=False),
                     f"[inference] SETTLES={'YES' if stats['settles'] else 'NO'} under the preregistered rule.",
                     '[observation] step elai_step likelihood reported_increment absolute_increment relative_increment'])
        for t in config['trace']:
            rows.append(' '.join(str(t[k]) for k in ('step', 'elai_step', 'likelihood', 'reported_increment',
                                                     'absolute_increment', 'relative_increment')))
        rows.extend('[observation] ERROR: ' + e for e in config['errors'])
    return '\n'.join(rows) + '\n'


def record_provenance(status_path: Path, output: Path, results: list[dict]) -> None:
    sys.path.insert(0, str(ROOT))
    from scripts import provenance
    # Operator forbids ALL git commands, including read-only provenance probes.
    provenance.get_git_info = lambda: {'disabled': 'operator prohibits all git commands'}
    recorder = provenance.RunRecorder(output.stem, command=sys.argv)
    for path in (Path(__file__), ROOT / 'scripts/elai_settling_sweep.sh',
                 ROOT / 'docs/elai-settling-sweep.md', BINARY,
                 BASE / 'work_long/panel.inp', BASE / 'work_long/panel.pos', status_path):
        recorder.record_input(path)
    recorder.set_seed(SEED)
    recorder.record_parameter('configurations', [{k: v for k, v in c.items() if k != 'trace'} for c in results])
    recorder.record_parameter('decision_rule', dict(decay_ratio=DECAY_RATIO, negative_limit=NEGATIVE_LIMIT, tolerance=TOLERANCE))
    recorder.record_output(output)
    for c in results:
        recorder.record_output(c['log'])
        work = BASE / (Path(c['log']).name.removesuffix('.stdout.txt') + '.work')
        for path in sorted((work / 'output').glob('*')):
            if path.is_file() and not path.is_symlink():
                recorder.record_output(path)
    recorder.save()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', action='store_true', help='print the runner configuration TSV')
    parser.add_argument('--status', type=Path, help='runner TSV: name steps returncode elapsed_seconds')
    parser.add_argument('--prefix', choices=('sweep', 'sweep_smoke'), default='sweep')
    parser.add_argument('--output', type=Path, help='new JSON report; refuses overwrite')
    parser.add_argument('--provenance', action='store_true')
    args = parser.parse_args()
    if args.plan:
        for config in CONFIGS:
            print('\t'.join(map(str, config + (SEED,))))
        return 0
    if args.status is None or args.output is None:
        parser.error('--status and --output are required unless --plan')
    if args.output.exists():
        parser.error('report exists; refusing overwrite')
    configs = {c[0]: c for c in CONFIGS}
    results, names = [], set()
    for line in args.status.read_text().splitlines():
        name, steps, returncode, elapsed = line.split('\t')
        if name not in configs or name in names:
            parser.error('unknown or repeated configuration')
        names.add(name)
        steps, returncode, elapsed = int(steps), int(returncode), int(elapsed)
        if steps != (5 if args.prefix == 'sweep_smoke' else configs[name][5]):
            parser.error('step count disagrees with registered plan')
        log = BASE / f'{args.prefix}_{name}.stdout.txt'
        text = log.read_text()
        errors = []
        try:
            trace = parse_trace(text)
        except ValueError as e:
            trace = []
            errors.append(str(e))
        if returncode not in (0, 1):
            errors.append(f'ELAI failed or timed out: returncode={returncode}')
        if '-elai: finished, for details see log:' not in text:
            errors.append('missing ELAI completion sentinel')
        if len(trace) != steps:
            errors.append(f'expected {steps} exact updates, found {len(trace)}')
        results.append(dict(name=name, configuration=dict(zip(('S', 'K', 'w', 'mg', 'planned_steps'), configs[name][1:])),
                            seed=SEED, log=str(log.relative_to(ROOT)), requested_steps=steps,
                            returncode=returncode, elapsed_seconds=elapsed, errors=errors, trace=trace,
                            statistics=decision(trace, not errors)))
    if not results:
        parser.error('empty status file')
    result = dict(validation_smoke=args.prefix == 'sweep_smoke', biological_use_authorized=False,
                  rule=dict(decay_ratio=DECAY_RATIO, negative_fraction_strictly_below=NEGATIVE_LIMIT,
                            relative_tolerance=TOLERANCE), configurations=results)
    with args.output.open('x') as f:
        json.dump(result, f, indent=2, allow_nan=False)
        f.write('\n')
    print(render(result), end='')
    if args.provenance:
        record_provenance(args.status, args.output, results)
    return int(any(c['errors'] for c in results))


if __name__ == '__main__':
    raise SystemExit(main())
