"""Software and full-panel artifact regressions for haplotype_background.py."""
import json
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from analysis import haplotype_background as h


class FakeGenotypes:
    def __init__(self, matrix):
        self.matrix = matrix
        self.n_samples = matrix.shape[1]

    def read_variants(self, indices):
        assert np.asarray(indices).dtype == np.uint32
        return self.matrix[np.asarray(indices)]


class BackgroundSoftwareTests(unittest.TestCase):
    def test_tested_variant_and_inclusive_radius_exclusion_invariance(self):
        centre = (h.LOCI['LCT'][1]+h.LOCI['LCT'][2])//2
        _, outerlo, outerhi, lo, hi = h.bounds('LCT', 200000)
        positions = np.array([lo-1, lo, centre, hi, hi+1, outerlo, outerhi])
        panel = pd.DataFrame({'index': np.arange(7), 'chrom': 2, 'pos': positions,
                              'ID': ['test_'+str(x) for x in positions], 'ref': 'C', 'alt': 'T'})
        rng = np.random.default_rng(3)
        g = rng.integers(0, 3, (7, 80)).astype(float)
        forbidden = h.excluded_positions(positions, lo, hi)
        self.assertEqual(forbidden.tolist(), [False, True, True, True, False, False, False])
        h.OUT.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=h.OUT) as scratch, patch.object(h, 'OUT', Path(scratch)):
            original = h.build_sketch(FakeGenotypes(g), panel, 'LCT')
            perturbed = g.copy(); perturbed[forbidden] = 2-g[forbidden]
            changed = h.build_sketch(FakeGenotypes(perturbed), panel, 'LCT')
            for width in (200000, 1000000):
                np.testing.assert_array_equal(original[0][width], changed[0][width])
                self.assertEqual(original[2][width], [0, 4, 5, 6])
            self.assertFalse(np.array_equal(original[3], changed[3]))

    def test_annotation_complete_exclusions_include_last_amy_copy_and_controls(self):
        h.initialize_intervals()
        self.assertEqual(h.EXCLUSIONS['AMY'], (1, 103915000, 104327000))
        for name, audit in h.EXCLUSION_AUDIT.items():
            self.assertEqual(audit['uncovered_core_transcripts'], 0)
            self.assertTrue(audit['assertions_passed'])
            t = pd.read_csv(h.OUT/f'{name}_core_transcripts.tsv', sep='\t')
            _, lo, hi = h.EXCLUSIONS[name]
            self.assertTrue((t.start.ge(lo) & t.end.le(hi)).all())
        amy = pd.read_csv(h.OUT/'AMY_core_transcripts.tsv', sep='\t')
        self.assertEqual(amy.end.max(), 104301314)
        self.assertGreater(h.EXCLUSION_AUDIT['AMY']['original_focal_transcripts_uncovered'], 0)

    def test_mean_imputation_and_invalid_dosage_are_not_observed_calls(self):
        g = np.array([[0, 1, 2, -9, np.nan, 3], [1, 1, 1, 1, 1, 1]], float)
        x, eligible, called = h.encode_block(g, np.array([0, 1], np.uint32))
        self.assertEqual(eligible.tolist(), [True, False])
        self.assertEqual(called[0].tolist(), [True, True, True, False, False, False])
        np.testing.assert_array_equal(x[3:], np.zeros((3, h.DIM)))

    def test_k_selection_does_not_use_outcome_and_ties_choose_smallest(self):
        self.assertEqual(h.choose_k({2: .45, 4: .46, 6: .2, 8: .1}), 2)
        self.assertEqual(h.choose_k({2: .1, 4: .4, 6: .3, 8: .2}), 4)

    def test_projection_is_batch_and_window_independent(self):
        ids = np.array([1, 5, 100], np.uint32)
        together = h.projection_weights(ids)
        np.testing.assert_array_equal(together, np.vstack([h.projection_weights([i]) for i in ids]))

    def test_temporal_floor_units_and_common_conditional_slope(self):
        rng = np.random.default_rng(7)
        n = 300
        labels = np.repeat([0, 1], n//2)
        time = np.tile(np.linspace(-8, -1, n//2), 2)
        y = .5+.025*time+.05*labels+rng.normal(0, .06, n)
        fitted = h.temporal_fit(y, time, labels)
        self.assertTrue(fitted['HC3_identifiable'])
        self.assertAlmostEqual(fitted['beta_per_kyr'], .025, delta=.005)
        self.assertLess(fitted['MDE80_HC3_per_kyr'], .025)
        self.assertEqual(fitted['MDE80_HC3_2kyr'], 2*fitted['MDE80_HC3_per_kyr'])
        expected = 2.802*np.sqrt(.25/(2*sum((time[:150]-time[:150].mean())**2)*2))
        self.assertAlmostEqual(h.design_floor(time, labels, .5)['MDE80_ideal_per_kyr'], expected)
        self.assertFalse(h.design_floor(time[:10], labels[:10], .5)['identifiable'])

    def test_negative_control_reproducible(self):
        labels = np.repeat(np.arange(4), 100)
        self.assertEqual(h.permutation_control(labels), h.permutation_control(labels))


class BackgroundArtifactTests(unittest.TestCase):
    def test_corrected_run_reproduces_scientific_results_and_retains_original_prereg(self):
        import hashlib
        original = h.OUT/'prereg.original.md'
        initial = json.loads((h.OUT/'invalid_exclusion_run/summary.json').read_text())
        self.assertEqual(hashlib.sha256(original.read_bytes()).hexdigest(), initial['prereg_sha256'])
        for name in h.LOCI:
            first = json.loads((h.OUT/f'{name}_summary.corrected_first.json').read_text())
            final = json.loads((h.OUT/f'{name}_summary.json').read_text())
            for r in [first, final]:
                r.pop('elapsed_seconds'); r.pop('peak_RSS_MiB')
            self.assertEqual(first, final)

    def test_full_panel_prereg_exclusion_gates_and_fail_closed(self):
        path = h.OUT/'summary.json'
        self.assertTrue(path.exists(), 'Run biological full-panel command before artifact validation')
        s = json.loads(path.read_text())
        self.assertTrue(s['status'].startswith('complete'))
        self.assertEqual(s['panel']['samples'], 16339)
        self.assertEqual(s['panel']['nominal_chromosome_slots'], 32678)
        self.assertEqual(s['panel']['variants'], 52382872)
        self.assertLess(s['prereg_mtime_ns'], datetime.fromisoformat(s['started_utc']).timestamp()*1e9)
        self.assertEqual(set(s['loci']), set(h.LOCI))
        for name, r in s['loci'].items():
            assignments = pd.read_csv(h.OUT/f'{name}_assignments.tsv', sep='\t')
            np.testing.assert_array_equal(assignments['index'], np.arange(16339))
            self.assertFalse(assignments.isna().any().any())
            for w in (200000, 1000000):
                features = pd.read_csv(h.OUT/f'{name}_{w}_features.tsv', sep='\t')
                self.assertFalse(features.excluded.any())
                _, _, _, lo, hi = h.bounds(name, w)
                self.assertFalse(h.excluded_positions(features.pos, lo, hi).any())
                self.assertEqual(len(features), r['coverage'][str(w)]['features'])
            self.assertEqual(r['interval_audit']['uncovered_core_transcripts'], 0)
            self.assertTrue(r['interval_audit']['assertions_passed'])
            self.assertEqual(set(r['gates']), set('abcd'))
            self.assertEqual(r['all_gates_pass'], all(r['gates'].values()))
            if not r['all_gates_pass']:
                self.assertTrue(r['outcomes']['status'].startswith('STOPPED'))
                self.assertFalse((h.OUT/f'{name}_temporal.tsv').exists())
        self.assertIn('outcome_blind_design_floor', s['loci']['AMY'])


if __name__ == '__main__': unittest.main()
