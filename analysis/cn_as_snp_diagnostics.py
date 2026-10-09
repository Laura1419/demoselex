#!/usr/bin/env python3
"""Labelled secondary diagnostics; never replace the primary CN-as-SNP fits."""
import argparse
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from cn_as_snp import (OUT, SEED, fit_binomial, fit_row, empirical_p, save_json,
                       AkbariGeno, parse_psam, paths, iter_pvar, AkbariGLMM)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--h9-max-iter', type=int, default=400)
    ap.add_argument('--permutations', type=int, default=199)
    args = ap.parse_args()
    primary = json.loads((OUT / 'summary.json').read_text())
    cohort = pd.read_csv(OUT / 'primary_cohort.tsv', sep='\t')
    synthetic = pd.read_csv(OUT / 'synthetic_loci.tsv', sep='\t')
    assert list(cohort['sample']) == list(synthetic['sample'])
    y = synthetic['H9'].to_numpy(float)
    t = cohort.age.to_numpy(float)
    data = pd.DataFrame({'time_bp': t, 'alt_count': y, 'total_count': 2})
    k = np.zeros((len(t), len(t)))
    diagnostic = {'kind': 'secondary_not_replacement', 'h9_max_iter': args.h9_max_iter}
    try:
        r = AkbariGLMM().fit(data, K=k, max_iter=args.h9_max_iter)
        r['X_equivalent'] = r['z_score'] / np.sqrt(2.28)
        rng = np.random.default_rng(SEED)
        strata = [np.flatnonzero(cohort.p2.to_numpy() == group) for group in cohort.p2.unique()]
        null, failed = [], 0
        with (OUT / 'H9_GRM_off_diagnostic_permutations.tsv').open('w') as f:
            f.write('permutation\tz_score\terror\n')
            for b in range(args.permutations):
                p = y.copy()
                for ix in strata:
                    p[ix] = y[rng.permutation(ix)]
                permdata = data.copy()
                permdata['alt_count'] = p
                try:
                    z = AkbariGLMM().fit(permdata, K=k, max_iter=args.h9_max_iter)['z_score']
                    null.append(z)
                    f.write(f'{b}\t{z}\t\n')
                except Exception as e:
                    failed += 1
                    f.write(f'{b}\t\t{e}\n')
                f.flush()
        r.update(empirical_P=empirical_p(r['z_score'], null),
                 permutations_requested=args.permutations, permutations_completed=len(null),
                 permutations_failed=failed,
                 null_abs_Z_quantiles=np.quantile(abs(np.asarray(null)), [.5, .9, .95, .99]).tolist() if null else [])
        diagnostic['H9_GRM_off'] = r
    except Exception as e:
        diagnostic['H9_GRM_off'] = {'converged': False, 'error': str(e)}
    save_json('diagnostics.json', diagnostic)
    print('H9 diagnostic', diagnostic['H9_GRM_off'], flush=True)

    ids = parse_psam(paths()[2])
    lookup = {s: i for i, s in enumerate(ids)}
    table = pd.read_excel('data/raw/akbari/41586_2026_10358_MOESM3_ESM.xlsx', sheet_name=0, header=1)
    datecol = 'Date mean in BP in years before 1950 CE'
    table['age'] = pd.to_numeric(table[datecol], errors='coerce')
    table['sample'] = table['Genetic version identifier'].astype(str)
    eligible = table[table['sample'].isin(lookup) & table.age.ge(0) &
                     table['Data source'].astype(str).str.contains('ancient', case=False) &
                     table['In unrelated set?'].eq(1)].drop_duplicates('sample').copy()
    rng = np.random.default_rng(SEED)
    secondary = eligible.iloc[rng.choice(len(eligible), min(600, len(eligible)), replace=False)].copy()
    secondary['pgen_full_index'] = secondary['sample'].map(lookup)
    secondary = secondary.sort_values('pgen_full_index').reset_index(drop=True)
    secondary[['sample', 'age', 'pgen_full_index']].to_csv(OUT / 'secondary_control_cohort.tsv', sep='\t', index=False)
    diagnostic['secondary_cohort'] = {'eligible': len(eligible), 'selected': len(secondary),
                                       'selection': 'seeded genotype-blind unrelated ancient Table S1 exact-ID sample'}
    needed = {primary['references'][n]['ID'] for n in ['LCT', 'SLC45A2']}
    positions = {}
    for i, (_, _, vid, _, _) in enumerate(iter_pvar(paths()[1])):
        if vid in needed:
            positions[vid] = i
            if len(positions) == len(needed):
                break
    bg = pd.read_csv(OUT / 'background_markers.tsv', sep='\t')
    with AkbariGeno(sample_subset=secondary.pgen_full_index.to_numpy(dtype=np.uint32)) as geno:
        g = geno.read_variants(bg['index'].to_numpy(dtype=np.uint32)).T.astype(float)
        g[(g < 0) | (g > 2)] = np.nan
        rows = []
        for name in ['LCT', 'SLC45A2']:
            ref = primary['references'][name]
            chrom = int(ref['CHROM'])
            k = AkbariGLMM.compute_grm(g, marker_chromosomes=bg.chrom.to_numpy(), target_chromosome=chrom)
            y = geno.read_hardcalls(np.array([positions[ref['ID']]], dtype=np.uint32))[0].astype(float)
            y[y < 0] = np.nan
            for model, mat in [('GRM_on', k), ('GRM_off', np.zeros_like(k))]:
                r = fit_row(name, model, y, secondary.age.to_numpy(float), mat)
                r['released_X'] = ref['X']
                rows.append(r)
                print('secondary control', r, flush=True)
        diagnostic['secondary_controls'] = rows
    save_json('diagnostics.json', diagnostic)


if __name__ == '__main__':
    main()
