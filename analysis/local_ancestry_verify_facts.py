#!/usr/bin/env python3
"""Recheck prior/global AMY facts without testing outcomes against failed local calls."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.local_ancestry import OUT, fit_frequency, save_json
from analysis.load_aadr import load_bundle
from analysis.load_akbari import iter_chunks


def main():
    meta = pd.read_csv(OUT/'scratch/ancient_metadata_for_supplement.tsv', sep='\t')
    cn = pd.read_csv('data/raw/amy/copy_number_table_for_supplement.tsv', sep='\t')
    ancient = cn[cn.type == 'Ancient'].merge(meta[['sample', 'age']], on='sample', validate='one_to_one')
    ancient = ancient[np.isfinite(ancient.age) & np.isfinite(ancient.AMY1)]
    means = ancient.groupby('p2').agg(n=('sample', 'size'), mean_AMY1=('AMY1', 'mean')).reset_index()
    means.to_csv(OUT/'verified_CN_means.tsv', sep='\t', index=False)
    _, inds, anno = load_bundle()
    by_key = {}
    for ind in inds:
        if ind.id not in anno:
            continue
        a = anno[ind.id]
        for key in {ind.id, ind.id.split('.')[0], a.individual_id}:
            by_key.setdefault(key, []).append(a)
    mapped = 0
    mapping = []
    for r in ancient.itertuples():
        aa = by_key.get(r.sample, [])
        bios = {a.individual_id for a in aa}
        ok = len(bios) == 1 and any(a.has_date and a.has_coords for a in aa)
        mapped += ok
        mapping.append(dict(sample=r.sample, unique_biological_match=len(bios) == 1, dated_coordinates=ok))
    pd.DataFrame(mapping).to_csv(OUT/'verified_CN_mapping.tsv', sep='\t', index=False)
    previous = Path('data/derived/amy_binarized')
    cohort = pd.read_csv(previous/'primary_cohort.tsv', sep='\t')
    dos = pd.read_csv(previous/'dosages.tsv', sep='\t')
    dd = cohort.merge(dos, on='sample', validate='one_to_one')
    assert len(dd) == 225
    old = json.loads((previous/'summary.json').read_text())
    baseline = {}
    for cut in ['primary_le1', 'sensitivity_le2']:
        y = dd[cut].to_numpy()/2
        time = -dd.age.to_numpy()/1000
        cov = pd.get_dummies(dd.p2, drop_first=True).to_numpy(float)
        fit = fit_frequency(y, time, cov)
        mu = pd.Series(y).groupby(dd.p2).transform('mean').to_numpy()
        centered = time-time.mean()
        beta_between = float(centered@mu/(centered@centered))
        beta_pooled = float(centered@y/(centered@centered))
        time_within = time-pd.Series(time).groupby(dd.p2).transform('mean').to_numpy()
        within_contrib = fit['beta_time']*float(time_within@time_within/(centered@centered))
        carriage = abs(beta_between)/(abs(beta_between)+abs(within_contrib))
        np.testing.assert_allclose(fit['beta_time'], old['channels'][cut]['beta_within'], atol=1e-10)
        np.testing.assert_allclose(fit['MDE80_HC3'], old['channels'][cut]['MDE80_HC3_per_kyr'], atol=1e-10)
        np.testing.assert_allclose(carriage, old['channels'][cut]['carriage_absolute_share'], atol=1e-10)
        baseline[cut] = dict(**fit, beta_pooled=beta_pooled, beta_between=beta_between,
                             within_pooled_contribution=within_contrib, carriage_share=carriage)
    results = old['results']
    primary = {r['model']: r for r in results if r['locus'] == 'primary_le1'}
    le2 = [r for r in results if r['locus'] == 'sensitivity_le2' and r['model'] == 'GRM_on'][0]
    # These two quantities are artifact checks, not a fresh GLMM fit.
    baseline['saved_GLMM'] = dict(primary_GRM_retention=abs(primary['GRM_on']['z_score']/primary['GRM_off']['z_score']),
                                   sensitivity_le2_X_equivalent=abs(le2['X_equivalent']))
    save_json('verified_global_baseline.json', baseline)
    windows = {'AMY': (1, 103500000, 104500000), 'LCT': (2, 136300000, 136900000),
               'FADS': (11, 61300000, 61700000), 'SLC45A2': (5, 33800000, 34200000),
               'SLC24A5': (15, 48200000, 48600000)}
    audit = {k: dict(n=0, max_abs_X=0., above_5_45=0) for k in windows}
    for chunk in iter_chunks(chunksize=500000, usecols=['CHROM', 'POS', 'X']):
        for name, (c, lo, hi) in windows.items():
            d = chunk[chunk.CHROM.eq(c) & chunk.POS.between(lo, hi)]
            if len(d):
                audit[name]['n'] += len(d)
                audit[name]['max_abs_X'] = max(audit[name]['max_abs_X'], float(d.X.abs().max()))
                audit[name]['above_5_45'] += int(d.X.abs().gt(5.45).sum())
    save_json('verified_akbari_windows.json', audit)
    save_json('verified_fact_counts.json', dict(dated_CN_ancients=len(ancient), CN_AADR_dated_coords=mapped,
                                                previous_labellable=old['cohort']['labellable'], previous_pgen_matched=len(dd)))
    print('Raw CN means and previous conditional/channel coefficients reproduced.', flush=True)
    print('CN AADR dated-coordinate matches:', mapped, 'of', len(ancient), flush=True)
    print('Raw released Akbari audit:', audit, flush=True)
    from scripts import provenance as p
    p.get_git_info = lambda: {'omitted': 'operator prohibits all git commands'}
    rec = p.RunRecorder('local-ancestry-carried-facts')
    for path in [Path(__file__), Path('data/raw/amy/copy_number_table_for_supplement.tsv'),
                 previous/'summary.json', previous/'primary_cohort.tsv', previous/'dosages.tsv']:
        rec.record_input(path)
    large = Path('data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz')
    rec.manifest['inputs'][str(large)] = dict(size=large.stat().st_size, mtime_ns=large.stat().st_mtime_ns)
    for path in OUT.glob('verified_*'):
        rec.record_output(path)
    rec.save()


if __name__ == '__main__':
    main()
