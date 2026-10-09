#!/usr/bin/env python3
"""Bounded, preregistered same-cohort CN-as-SNP comparison; validation first."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import zipfile

import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'analysis'))
sys.path.insert(0, str(ROOT / 'sim'))
from load_akbari_geno import AkbariGeno, parse_psam, paths, iter_pvar
from load_akbari import iter_chunks
from akbari_glmm import AkbariGLMM

OUT = Path('data/derived/cn_as_snp')
CF = 2.28
SEED = 20261009
HAPS = ['H1^a', 'H5', 'H7', 'H9', 'H3^r']
PREREG = Path('docs/cn-as-snp-prereg.md')


def ordering_gate(rows):
    """A failed/absent fit never passes the prespecified relative gate."""
    on = {r['locus']: r for r in rows if r['model'] == 'GRM_on'}
    required = ['LCT', 'SLC45A2', 'AMY_real_SNP']
    if any(k not in on or not on[k].get('converged', False) for k in required):
        return False
    weak = abs(on['AMY_real_SNP']['z_score'])
    return all(abs(on[k]['z_score']) >= 2 and abs(on[k]['z_score']) >= 1.25 * weak
               for k in ['LCT', 'SLC45A2'])


def completion_status(rows, requested):
    successful = [r for r in rows if r.get('converged')]
    if any(r.get('permutations_completed', 0) + r.get('permutations_failed', 0) < requested
           for r in successful):
        return 'partial_permutation_budget'
    if any(not r.get('converged') or r.get('permutations_failed', 0) for r in rows):
        return 'complete_with_fit_failures'
    return 'complete'


def empirical_p(observed, null):
    a = np.asarray(null, float)
    a = a[np.isfinite(a)]
    return float((1 + (abs(a) >= abs(observed)).sum()) / (1 + len(a))) if len(a) else None


def match_samples(ancient, ids):
    exact = {s: i for i, s in enumerate(ids)}
    bio = defaultdict(list)
    for i, s in enumerate(ids):
        bio[str(s).split('.')[0]].append(i)
    rows = []
    for _, row in ancient.iterrows():
        s = str(row['sample'])
        if s in exact:
            ix, method = exact[s], 'exact'
        elif len(bio[s.split('.')[0]]) == 1:
            ix, method = bio[s.split('.')[0]][0], 'unique_prefix'
        else:
            ix, method = None, 'ambiguous' if bio[s.split('.')[0]] else 'absent'
        r = row.to_dict()
        r.update(pgen_full_index=ix, pgen_id=ids[ix] if ix is not None else None,
                 match_method=method)
        rows.append(r)
    return pd.DataFrame(rows)


def extract_tables():
    archive = Path('data/raw/amy/amylase_diversity_project-v1.1.zip')
    tables = {}
    with zipfile.ZipFile(archive) as z:
        for suffix in ['haplotype_structure_for_supplement.tsv',
                       'ancient_metadata_for_supplement.tsv',
                       'copy_number_table_for_supplement.tsv']:
            names = [n for n in z.namelist() if n.endswith('/' + suffix)]
            if len(names) != 1:
                raise ValueError(f'ambiguous archive member: {suffix}')
            raw = z.read(names[0])
            (OUT / suffix).write_bytes(raw)
            tables[suffix] = pd.read_csv(OUT / suffix, sep='\t')
    h = tables['haplotype_structure_for_supplement.tsv']
    a = h[h.type.eq('Ancient')].copy()
    if a['sample'].duplicated().any():
        raise ValueError('duplicate ancient diplotypes')
    meta = tables['ancient_metadata_for_supplement.tsv']
    a = a.merge(meta[['sample', 'age', 'p2_l0']], on='sample', validate='one_to_one')
    a['age'] = pd.to_numeric(a.age, errors='coerce')
    a = a[a.age.notna() & a.H1.notna() & a.H2.notna()].copy()
    return a, tables['copy_number_table_for_supplement.tsv']


def references():
    canonical = []
    amy = []
    for chunk in iter_chunks(usecols=['ID', 'CHROM', 'POS', 'REF', 'ALT', 'X', 'AF']):
        mask = ((chunk.CHROM.eq(2) & chunk.POS.eq(136608646)) |
                (chunk.CHROM.eq(5) & chunk.POS.eq(33951693)))
        canonical.extend(chunk[mask].to_dict('records'))
        amy.extend(chunk[chunk.CHROM.eq(1) & chunk.POS.between(103500000, 104500000)].to_dict('records'))
    selected = {}
    for label, chrom in [('LCT', 2), ('SLC45A2', 5)]:
        rows = [r for r in canonical if int(r['CHROM']) == chrom]
        if len(rows) != 1:
            raise ValueError(f'exact canonical SNP not unique: {label} {rows}')
        selected[label] = rows[0]
    if not amy:
        raise ValueError('no released AMY SNPs')
    selected['AMY_real_SNP'] = max(amy, key=lambda r: abs(r['X']))
    window = {'n_variants': len(amy), 'max_abs_X': max(abs(r['X']) for r in amy),
              'n_above_5_45': sum(abs(r['X']) > 5.45 for r in amy)}
    return selected, window


def map_markers(n_variants, selected):
    grid = set(np.linspace(0, n_variants - 1, min(10000, n_variants), dtype=np.uint32).tolist())
    needed = {v['ID'] for v in selected.values()}
    background = []
    positions = {}
    for i, (chrom, pos, vid, ref, alt) in enumerate(iter_pvar(paths()[1])):
        if vid in needed:
            positions[vid] = i
        if i in grid and chrom in {str(x) for x in range(1, 23)} and len(ref) == len(alt) == 1:
            background.append((i, int(chrom), pos, vid))
    if set(positions) != needed:
        raise ValueError('canonical or AMY real SNP missing from pvar')
    return background, positions


def fit_binomial(y, dates, k):
    keep = np.isfinite(y) & (y >= 0) & (y <= 2)
    data = pd.DataFrame({'time_bp': dates[keep], 'alt_count': y[keep], 'total_count': 2})
    r = AkbariGLMM().fit(data, K=k[np.ix_(keep, keep)])
    r['n'] = int(keep.sum())
    return r


def fit_continuous(y, dates, k):
    keep = np.isfinite(y)
    y, t, k = y[keep], dates[keep], k[np.ix_(keep, keep)]
    scale = t.std()
    x = np.column_stack([np.ones(len(y)), (t - t.mean()) / scale])
    var = np.var(y)
    if len(y) < 4 or var <= 0 or scale <= 0:
        raise ValueError('insufficient quantitative variation')
    def objective(v):
        return AkbariGLMM._gls(y, x, v[0] * k + v[1] * np.eye(len(y)))[0]
    opt = minimize(objective, [var * .2, var * .8], method='L-BFGS-B',
                   bounds=[(0, None), (1e-8, None)])
    if not opt.success:
        raise RuntimeError(f'Gaussian REML failed: {opt.message}')
    _, beta, cov, _ = AkbariGLMM._gls(y, x, opt.x[0] * k + opt.x[1] * np.eye(len(y)))
    return {'z_score': float(beta[1] / np.sqrt(cov[1, 1])),
            'beta_time': float(beta[1] / scale), 'se_beta_time': float(np.sqrt(cov[1, 1]) / scale),
            'tau1': float(opt.x[0]), 'residual_variance': float(opt.x[1]),
            'n': int(len(y)), 'converged': True}


def fit_row(locus, model, y, dates, k, continuous=False):
    try:
        r = (fit_continuous if continuous else fit_binomial)(y, dates, k)
        r.update(X_equivalent=r['z_score'] / np.sqrt(CF),
                 clears_3_16=bool(abs(r['z_score']) / np.sqrt(CF) > 3.16),
                 clears_5_45=bool(abs(r['z_score']) / np.sqrt(CF) > 5.45))
    except Exception as e:
        r = {'converged': False, 'error': str(e)}
    return dict(locus=locus, model=model, **r)


def save_json(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, default=str) + '\n')


def record_manifest(args, start):
    # The repository harness executes git in its constructor. Override ONLY that
    # metadata provider to respect the explicit no-git task boundary.
    sys.path.insert(0, str(ROOT))
    from scripts import provenance
    provenance.get_git_info = lambda: {'commit': 'not inspected: no-git directive',
                                       'branch': 'not inspected', 'dirty': 'unknown'}
    recorder = provenance.RunRecorder('cn-as-snp', command=[sys.executable] + sys.argv)
    recorder.record_parameter('configuration', vars(args))
    recorder.record_parameter('started_utc', start)
    recorder.record_parameter('prereg_sha256', hashlib.sha256(PREREG.read_bytes()).hexdigest())
    recorder.set_seed(SEED)
    for p in [PREREG, Path('analysis/cn_as_snp.py'), Path('sim/akbari_glmm.py'), paths()[2],
              Path('data/raw/amy/amylase_diversity_project-v1.1.zip')]:
        recorder.record_input(p)
    # Large immutable genotype/stat files identified by size/mtime; no falsely
    # asserted content hash and no lengthy rehash of the 25GB input.
    for p in [paths()[0], paths()[1], Path('data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz')]:
        st = p.stat()
        recorder.record_parameter(str(p), {'size': st.st_size, 'mtime_ns': st.st_mtime_ns,
                                          'sha256': 'not computed within bounded run'})
    for p in OUT.glob('*'):
        if p.is_file():
            recorder.record_output(p)
    recorder.save()


def run(args):
    if not PREREG.exists():
        raise RuntimeError('preregistration required before fitting')
    OUT.mkdir(parents=True, exist_ok=True)
    start = datetime.now(timezone.utc).isoformat()
    deadline = time.monotonic() + args.seconds
    state = {'started_utc': start, 'status': 'validation_in_progress', 'seed': SEED,
             'prereg_mtime_ns': PREREG.stat().st_mtime_ns, 'results': []}
    save_json('summary.json', state)
    ancient, copies = extract_tables()
    ids = parse_psam(paths()[2])
    mapped = match_samples(ancient, ids)
    mapped.to_csv(OUT / 'sample_mapping.tsv', sep='\t', index=False)
    cohort = mapped[mapped.pgen_full_index.notna()].copy()
    if cohort.pgen_full_index.duplicated().any():
        raise ValueError('multiple AMY individuals map to one genotype row')
    cohort = cohort.sort_values('pgen_full_index').reset_index(drop=True)
    cohort.to_csv(OUT / 'primary_cohort.tsv', sep='\t', index=False)
    state['cohort'] = {'labellable_dated_ancients': len(ancient), 'matched': len(cohort),
                       'match_rate': len(cohort) / len(ancient),
                       'mapping_methods': mapped.match_method.value_counts().to_dict()}
    print('cohort', state['cohort'], flush=True)
    selected, state['released_AMY_window'] = references()
    state['references'] = selected
    save_json('summary.json', state)
    with AkbariGeno(sample_subset=cohort.pgen_full_index.to_numpy(dtype=np.uint32)) as geno:
        bg, pos = map_markers(geno.n_variants, selected)
        bg = pd.DataFrame(bg, columns=['index', 'chrom', 'pos', 'ID'])
        g = geno.read_variants(bg['index'].to_numpy(dtype=np.uint32)).T.astype(float)
        g[(g < 0) | (g > 2)] = np.nan
        call = np.isfinite(g).mean(axis=0)
        p = np.nanmean(g, axis=0) / 2
        ok = (call >= .95) & (np.minimum(p, 1 - p) >= .05)
        g, bg = g[:, ok], bg[ok].reset_index(drop=True)
        bg.to_csv(OUT / 'background_markers.tsv', sep='\t', index=False)
        state['GRM'] = {'markers': len(bg), 'chromosome_counts': bg.chrom.value_counts().to_dict(),
                        'method': 'vanraden', 'minimum_call_rate': .95, 'minimum_MAF': .05,
                        'individual_level': True, 'LOCO': True}
        if len(cohort) < 40 or len(bg) < 100 or bg.chrom.nunique() < 10:
            state.update(status='stopped_insufficient_genomewide_GRM', gate_passed=False)
            save_json('summary.json', state)
            record_manifest(args, start)
            return state
        dates = cohort.age.to_numpy(float)
        zeros = np.zeros((len(cohort), len(cohort)))
        cache = {}
        def grm(chrom):
            if chrom not in cache:
                cache[chrom] = AkbariGLMM.compute_grm(g, marker_chromosomes=bg.chrom.to_numpy(),
                                                     target_chromosome=chrom)
                np.save(OUT / f'GRM_LOCO_chr{chrom}.npy', cache[chrom])
            return cache[chrom]
        phenotypes = {}
        chroms = {}
        for name, ref in selected.items():
            y = geno.read_hardcalls(np.array([pos[ref['ID']]], dtype=np.uint32))[0].astype(float)
            y[y < 0] = np.nan
            dose = geno.read_variants(np.array([pos[ref['ID']]], dtype=np.uint32))[0]
            state.setdefault('fractional_control_dosages', {})[name] = int(((dose >= 0) & (abs(dose - np.round(dose)) > .001)).sum())
            phenotypes[name], chroms[name] = y, int(ref['CHROM'])
            for model, k in [('GRM_on', grm(chroms[name])), ('GRM_off', zeros)]:
                r = fit_row(name, model, y, dates, k)
                r['released_X'] = float(ref['X'])
                state['results'].append(r)
                save_json('summary.json', state)
                print(json.dumps(r), flush=True)
        state['gate_passed'] = ordering_gate(state['results'])
        if not state['gate_passed']:
            state['status'] = 'stopped_ordering_gate_failed'
            save_json('summary.json', state)
            record_manifest(args, start)
            return state
        state['status'] = 'counterfactual_in_progress'
        save_json('summary.json', state)
        for h in HAPS:
            phenotypes[h] = (cohort.H1.eq(h).astype(int) + cohort.H2.eq(h).astype(int)).to_numpy(float)
            chroms[h] = 1
        cn = copies[copies.type.eq('Ancient')][['sample', 'AMY1']]
        cohort_cn = cohort[['sample']].merge(cn, how='left', validate='one_to_one')
        phenotypes['AMY1_continuous'] = cohort_cn.AMY1.to_numpy(float)
        chroms['AMY1_continuous'] = 1
        for c in [3, 4, 6]:
            marker = bg[bg.chrom.eq(c)].iloc[0]
            name = f'background_chr{c}'
            y = geno.read_hardcalls(np.array([marker['index']], dtype=np.uint32))[0].astype(float)
            y[y < 0] = np.nan
            phenotypes[name], chroms[name] = y, c
        descriptions = []
        bins = pd.cut(cohort.age, [0, 2000, 4000, 6000, 8000, 10000, np.inf], right=False)
        for h in HAPS:
            y = phenotypes[h]
            state.setdefault('dosage_counts', {})[h] = {str(i): int((y == i).sum()) for i in range(3)}
            for dimension, labels in [('ancestry', cohort.p2), ('date_bin', bins)]:
                for group in labels.dropna().unique():
                    ix = (labels == group).to_numpy()
                    descriptions.append({'haplotype': h, 'dimension': dimension, 'group': str(group),
                                         'n': int(ix.sum()), 'frequency': float(y[ix].sum() / (2 * ix.sum()))})
        pd.DataFrame(descriptions).to_csv(OUT / 'frequencies.tsv', sep='\t', index=False)
        pd.DataFrame({h: phenotypes[h] for h in HAPS} | {'sample': cohort['sample'],
                     'AMY1_continuous': phenotypes['AMY1_continuous']}).to_csv(OUT / 'synthetic_loci.tsv', sep='\t', index=False)
        for name, y in phenotypes.items():
            if name in selected:
                continue
            for model, k in [('GRM_on', grm(chroms[name])), ('GRM_off', zeros)]:
                r = fit_row(name, model, y, dates, k, name == 'AMY1_continuous')
                lct = next(x for x in state['results'] if x['locus'] == 'LCT' and x['model'] == model)
                if r.get('converged') and lct.get('converged'):
                    r['AMY_to_LCT_abs_Z_ratio'] = abs(r['z_score'] / lct['z_score'])
                    r['extrapolated_abs_X'] = abs(selected['LCT']['X']) * r['AMY_to_LCT_abs_Z_ratio']
                state['results'].append(r)
                save_json('summary.json', state)
                print(json.dumps(r), flush=True)
        rng = np.random.default_rng(SEED)
        strata = [np.flatnonzero(cohort.p2.to_numpy() == group) for group in cohort.p2.unique()]
        with (OUT / 'permutations.tsv').open('w') as f:
            f.write('locus\tmodel\tpermutation\tz_score\terror\n')
            # Round-robin allows every phenotype/model some null evidence before timeout.
            nulls = {(r['locus'], r['model']): [] for r in state['results'] if r.get('converged')}
            failed = {key: 0 for key in nulls}
            for b in range(args.permutations):
                if time.monotonic() >= deadline:
                    break
                for name, y in phenotypes.items():
                    perm = y.copy()
                    # Preserve missingness in place for fair control support.
                    for ix in strata:
                        ix = ix[np.isfinite(y[ix])]
                        perm[ix] = y[rng.permutation(ix)]
                    for model, k in [('GRM_on', grm(chroms[name])), ('GRM_off', zeros)]:
                        key = (name, model)
                        if key not in nulls or time.monotonic() >= deadline:
                            continue
                        r = fit_row(name, model, perm, dates, k, name == 'AMY1_continuous')
                        if r.get('converged'):
                            nulls[key].append(r['z_score'])
                        else:
                            failed[key] += 1
                        f.write(f"{name}\t{model}\t{b}\t{r.get('z_score', '')}\t{r.get('error', '')}\n")
                        f.flush()
                for r in state['results']:
                    key = (r['locus'], r['model'])
                    if key in nulls:
                        a = nulls[key]
                        r.update(permutations_requested=args.permutations, permutations_completed=len(a),
                                 permutations_failed=failed[key], empirical_P=empirical_p(r['z_score'], a),
                                 null_abs_Z_quantiles=np.quantile(abs(np.asarray(a)), [.5, .9, .95, .99]).tolist() if a else [])
                save_json('summary.json', state)
                print('permutation round', b + 1, flush=True)
        state['status'] = completion_status(state['results'], args.permutations)
    save_json('summary.json', state)
    pd.DataFrame(state['results']).to_csv(OUT / 'statistics.tsv', sep='\t', index=False)
    record_manifest(args, start)
    return state


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--seconds', type=int, default=1400)
    ap.add_argument('--permutations', type=int, default=199)
    args = ap.parse_args()
    run(args)


if __name__ == '__main__':
    main()
