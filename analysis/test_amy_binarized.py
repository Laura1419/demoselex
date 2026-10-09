"""Source-named regressions and saved-artifact audit for amy_binarized.py."""
import json
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from amy_binarized import hap_copies, dosage, channels, CUTS, OUT, PRE, empirical_p


class AmyBinarizedTests(unittest.TestCase):
    def test_hap_copies_separates_amy1_from_amy2(self):
        self.assertEqual(hap_copies('H2A2B2'),2)
        self.assertEqual(hap_copies('H3^r'),3)
        self.assertEqual(hap_copies('H1^a'),1)
        for label in ['unknown',None,'H3unexpected']:
            with self.assertRaises(ValueError): hap_copies(label)

    def test_dosage_all_three_cuts(self):
        d=pd.DataFrame({'H1':['H1^a','H2A0','H3^r','H5'],'H2':['H1^a','H1^a','H2A2B2','H7']})
        np.testing.assert_array_equal(dosage(d,1),[0,1,2,2])
        np.testing.assert_array_equal(dosage(d,2),[0,0,1,2])
        np.testing.assert_array_equal(dosage(d,3),[0,0,0,2])

    def test_channels_sum_pooled_and_baseline_invariance(self):
        age=np.array([10000,9000,8000,7000,4000,3000,2000,1000])
        g=np.repeat(['old','young'],4); y=np.array([0,0,1,1,1,1,2,2])
        c=channels(age,y,g)
        self.assertAlmostEqual(c['beta_pooled'],c['beta_between']+c['within_pooled_contribution'])
        self.assertAlmostEqual(c['beta_total'],c['beta_between']+c['beta_within'])
        shifted=y+np.repeat([.2,-.2],4)
        self.assertAlmostEqual(channels(age,shifted,g)['beta_within'],c['beta_within'])
        self.assertGreater(c['MDE80_HC3_per_kyr'],0)

    def test_actual_dosages_and_frequencies(self):
        d=pd.read_csv(OUT/'primary_cohort.tsv',sep='\t'); y=pd.read_csv(OUT/'dosages.tsv',sep='\t')
        self.assertEqual(len(d),225); self.assertEqual(list(d['sample']),list(y['sample']))
        for name,cut in CUTS.items(): np.testing.assert_array_equal(y[name],dosage(d,cut))
        old=pd.read_csv('data/derived/cn_as_snp/synthetic_loci.tsv',sep='\t')
        np.testing.assert_array_equal(y.primary_le1,2-old['H1^a'])
        f=pd.read_csv(OUT/'frequencies.tsv',sep='\t')
        for _,r in f.iterrows():
            self.assertEqual(r.n,r.n0+r.n1+r.n2)
            self.assertEqual(r.B_copies,r.n1+2*r.n2)
            self.assertEqual(r.alleles,2*r.n)
            self.assertAlmostEqual(r.B_frequency,r.B_copies/r.alleles)
        for _,p in f.groupby(['scope','cut','dimension']):
            self.assertEqual(p.n.sum(),225 if p.scope.iloc[0]=='matched_225' else 288)

    def test_fits_and_prereg_timestamp(self):
        s=json.loads((OUT/'summary.json').read_text())
        self.assertLess(s['prereg_mtime_ns'],(OUT/'dosages.tsv').stat().st_mtime_ns)
        self.assertIn('MODELING CHOICE',PRE.read_text())
        self.assertTrue(s['ordering_gate_passed'])
        self.assertEqual(s['GRM_verified']['AMY_LOCO_markers'],1058)
        self.assertEqual(s['GRM_verified']['background_markers'],1133)
        old=json.loads(Path('data/derived/cn_as_snp/summary.json').read_text())
        for r in s['controls']:
            prior=next(q for q in old['results'] if q['locus']==r['locus'] and q['model']==r['model'])
            self.assertAlmostEqual(r['z_score'],prior['z_score'],places=8)
        for r in s['results']:
            self.assertTrue(r['converged']); self.assertEqual(r['n'],225)
            self.assertAlmostEqual(r['X_equivalent'],r['z_score']/np.sqrt(2.28))
        for model in ['GRM_on','GRM_off']:
            primary=next(r for r in s['results'] if r['locus']=='primary_le1' and r['model']==model)
            ancestral=next(r for r in old['results'] if r['locus']=='H1^a' and r['model']==model)
            self.assertAlmostEqual(primary['z_score'],-ancestral['z_score'],places=6)
        for c in s['channels'].values():
            self.assertAlmostEqual(c['beta_total'],c['beta_within']+c['beta_between'])
            self.assertAlmostEqual(c['beta_pooled'],c['within_pooled_contribution']+c['beta_between'])

    def test_permutation_artifacts(self):
        s=json.loads((OUT/'summary.json').read_text()); p=pd.read_csv(OUT/'permutations.tsv',sep='\t')
        obs={(r['locus'],r['model']):r['z_score'] for r in s['results']}
        self.assertEqual(len(s['nulls']),12)
        for r in s['nulls']:
            rows=p[p['null'].eq(r['null'])&p['cut'].eq(r['cut'])&p.model.eq(r['model'])]
            self.assertEqual(len(rows),r['successful']+r['failed'])
            vals=rows.loc[rows.converged,'z_score'].to_numpy()
            self.assertEqual(len(vals),r['successful'])
            self.assertAlmostEqual(empirical_p(obs[(r['cut'],r['model'])],vals),r['empirical_P'])
            self.assertEqual(len(rows),199)

    def test_genomewide_distribution_blocks(self):
        d=json.loads((OUT/'distribution.json').read_text()); b=pd.read_csv(OUT/'genomewide_block_maxima.tsv',sep='\t')
        self.assertEqual(d['n_finite'],9739624)
        self.assertEqual(d['window_audit']['AMY']['n'],2204)
        self.assertEqual(d['window_audit']['AMY']['strong_n'],0)
        for name,audit in d['window_audit'].items():
            self.assertEqual(audit['max_abs_X'],d['thresholds'][name]['abs_X'])
        self.assertEqual(len(b),d['all_blocks'])
        self.assertEqual(int(b.max_abs_X.gt(5.45).sum()),d['strong_max_blocks'])
        self.assertEqual(int(b.max_abs_X.between(3.16,5.45).sum()),d['band_max_blocks'])
        self.assertEqual(d['band_any_blocks']-d['overlap_blocks'],d['band_only_any_blocks'])
        for r in d['thresholds'].values(): self.assertAlmostEqual(r['fraction_exceed'],r['n_exceed']/d['n_finite'])


if __name__=='__main__': unittest.main()
