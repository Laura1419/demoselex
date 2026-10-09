#!/usr/bin/env python3
"""Prespecified secondary unrestricted-label null; primary stratified null retained."""
import json
import time
import numpy as np
import pandas as pd
from cn_as_snp import (OUT, SEED, HAPS, AkbariGeno, AkbariGLMM, paths, iter_pvar,
                       fit_row, empirical_p, save_json)


def main():
    # Fixed first-run snapshot prevents reading a half-written ongoing summary.
    primary = json.loads((OUT / 'first_run_summary.json').read_text())
    diag = json.loads((OUT / 'diagnostics.json').read_text())
    cohort = pd.read_csv(OUT / 'primary_cohort.tsv', sep='\t')
    syn = pd.read_csv(OUT / 'synthetic_loci.tsv', sep='\t')
    assert list(cohort['sample']) == list(syn['sample'])
    dates = cohort.age.to_numpy(float)
    phen = {h: syn[h].to_numpy(float) for h in HAPS + ['AMY1_continuous']}
    chroms = {h: 1 for h in phen}
    refs = primary['references']
    needed = {r['ID'] for r in refs.values()}
    positions = {}
    for i, (_, _, vid, _, _) in enumerate(iter_pvar(paths()[1])):
        if vid in needed:
            positions[vid] = i
            if len(positions) == len(needed):
                break
    bg = pd.read_csv(OUT / 'background_markers.tsv', sep='\t')
    with AkbariGeno(sample_subset=cohort.pgen_full_index.to_numpy(dtype=np.uint32)) as geno:
        for name, ref in refs.items():
            y = geno.read_hardcalls(np.array([positions[ref['ID']]], dtype=np.uint32))[0].astype(float)
            y[y < 0] = np.nan
            phen[name], chroms[name] = y, int(ref['CHROM'])
        for c in [3, 4, 6]:
            marker = bg[bg.chrom.eq(c)].iloc[0]
            y = geno.read_hardcalls(np.array([marker['index']], dtype=np.uint32))[0].astype(float)
            y[y < 0] = np.nan
            phen[f'background_chr{c}'], chroms[f'background_chr{c}'] = y, c
    covariance = {c: np.load(OUT / f'GRM_LOCO_chr{c}.npy') for c in set(chroms.values())}
    observed = {(r['locus'], r['model']): r for r in primary['results'] if r.get('converged')}
    observed[('H9', 'GRM_off')] = dict(locus='H9', model='GRM_off', **diag['H9_GRM_off'])
    nulls = {key: [] for key in observed}
    failed = {key: 0 for key in observed}
    rng = np.random.default_rng(SEED)
    deadline = time.monotonic() + 600
    with (OUT / 'unrestricted_permutations.tsv').open('w') as f:
        f.write('locus\tmodel\tpermutation\tz_score\terror\n')
        for b in range(199):
            if time.monotonic() > deadline:
                break
            for name, y in phen.items():
                ix = np.flatnonzero(np.isfinite(y))
                perm = y.copy()
                perm[ix] = y[rng.permutation(ix)]
                for model in ['GRM_on', 'GRM_off']:
                    key = (name, model)
                    if key not in observed or time.monotonic() > deadline:
                        continue
                    k = covariance[chroms[name]] if model == 'GRM_on' else np.zeros_like(covariance[chroms[name]])
                    if name == 'H9' and model == 'GRM_off':
                        try:
                            frame = pd.DataFrame({'time_bp': dates, 'alt_count': perm, 'total_count': 2})
                            r = AkbariGLMM().fit(frame, K=k, max_iter=400)
                        except Exception as e:
                            r = {'converged': False, 'error': str(e)}
                    else:
                        r = fit_row(name, model, perm, dates, k, name == 'AMY1_continuous')
                    if r.get('converged'):
                        nulls[key].append(r['z_score'])
                    else:
                        failed[key] += 1
                    f.write(f"{name}\t{model}\t{b}\t{r.get('z_score', '')}\t{r.get('error', '')}\n")
                    f.flush()
            rows = []
            for key, values in nulls.items():
                rows.append({'locus': key[0], 'model': key[1], 'empirical_P': empirical_p(observed[key]['z_score'], values),
                             'permutations_requested': 199, 'permutations_completed': len(values),
                             'permutations_failed': failed[key], 'max_iter': 400 if key == ('H9', 'GRM_off') else 40,
                             'null_abs_Z_quantiles': np.quantile(abs(np.asarray(values)), [.5, .9, .95, .99]).tolist() if values else []})
            save_json('unrestricted_null.json', {'null': 'secondary_unrestricted', 'results': rows})
            print('unrestricted round', b + 1, flush=True)


if __name__ == '__main__':
    main()
