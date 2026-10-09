"""Source-named software regressions for the Stage-1 ancestry approximation."""
import itertools
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from analysis.local_ancestry import (PAIRS, NAMES, COMPOSITION, source_frequency,
                                     transition, smooth_hmm, paint, fit_frequency)


class LocalAncestryTests(unittest.TestCase):
    def test_source_frequency_pseudohaploid_is_one_allele(self):
        g = np.array([[2, 3], [0, 2], [1, 1]])
        p, call, count = source_frequency(g, np.array([False, False, True]))
        np.testing.assert_allclose(p, [(2+.5)/(4+1), (2+.5)/(3+1)])
        np.testing.assert_allclose(call, [1, 2/3])
        np.testing.assert_array_equal(count, [3, 2])

    def test_transition_pair_alignment_and_row_normalization(self):
        tr = transition(4)
        np.testing.assert_allclose(tr.sum(axis=1), 1)
        aa = np.flatnonzero(np.all(PAIRS == [0, 0], axis=1))[0]
        ab = np.flatnonzero(np.all(PAIRS == [0, 1], axis=1))[0]
        bb = np.flatnonzero(np.all(PAIRS == [1, 1], axis=1))[0]
        self.assertAlmostEqual(tr[aa, ab]/tr[aa, aa], np.exp(-4))
        self.assertAlmostEqual(tr[aa, bb]/tr[aa, aa], np.exp(-8))
        self.assertTrue(np.allclose(COMPOSITION.sum(axis=1), 1))

    def test_hmm_forward_backward_and_viterbi_against_enumeration(self):
        rng = np.random.default_rng(16)
        logs = rng.normal(size=(1, 3, len(PAIRS)))
        post, path = smooth_hmm(logs, 2, viterbi=True)
        tr = transition(2)
        sequences = np.array(list(itertools.product(range(len(PAIRS)), repeat=3)))
        weights = np.array([np.exp(logs[0, np.arange(3), seq].sum())*tr[seq[0], seq[1]]*tr[seq[1], seq[2]]/len(PAIRS) for seq in sequences])
        expected = np.zeros((3, len(PAIRS)))
        for j in range(3):
            for k in range(len(PAIRS)):
                expected[j, k] = weights[sequences[:, j] == k].sum()/weights.sum()
        np.testing.assert_allclose(post[0], expected, atol=1e-12)
        np.testing.assert_array_equal(path[0], sequences[weights.argmax()])

    def test_missing_windows_are_flat_and_chromosomes_reset(self):
        panel = pd.DataFrame(dict(chrom=[1]*6+[2]*6, pos=list(range(0, 600000, 100000))*2))
        for k, s in enumerate(NAMES):
            panel['p_'+s] = .98 if k == 0 else .02
        g = np.array([[2]*6+[3]*6], dtype=np.uint8)
        fractions, support, windows, path = paint(g, np.array([False]), panel, 1000000, 8, True)
        self.assertGreater(fractions[0, 0, 0], .6)
        self.assertFalse(support[0, 1])
        # Equal unordered-pair prior has equal expected source fractions.
        np.testing.assert_allclose(fractions[0, 1], np.ones(5)/5)
        self.assertEqual(windows[1, 0], 2)

    def test_conditional_time_fit_and_detection_floor(self):
        rng = np.random.default_rng(100)
        time = np.linspace(-10, 0, 100)
        cov = rng.normal(size=(100, 1))
        y = .5+.02*time+.04*cov[:, 0]+rng.normal(0, .02, 100)
        result = fit_frequency(y, time, cov)
        self.assertTrue(result['identifiable'])
        self.assertAlmostEqual(result['beta_time'], .02, delta=.002)
        self.assertGreater(result['MDE80_HC3'], 0)
        self.assertLess(result['p'], .01)
        bad = fit_frequency(y, time, time[:, None])
        self.assertFalse(bad['identifiable'])

    def test_preregistered_artifacts_and_fail_closed(self):
        out = Path('data/derived/local_ancestry')
        if not (out/'validation.json').exists():
            self.skipTest('biological run not yet complete')
        import json
        validation = json.loads((out/'validation.json').read_text())
        self.assertLess(Path('docs/local-ancestry-prereg.md').stat().st_mtime_ns,
                        (out/'aim_panel.tsv').stat().st_mtime_ns)
        panel = pd.read_csv(out/'aim_panel.tsv', sep='\t')
        self.assertTrue((panel.divergence >= .45).all())
        self.assertTrue((panel.target_call_fraction >= .2).all())
        self.assertTrue((panel[[f'call_{s}' for s in NAMES]] >= .5).all().all())
        for _, group in panel.groupby('chrom'):
            self.assertTrue((np.diff(group.pos) >= 50000).all())
        self.assertEqual(validation['passed'], all(validation['gates'].values()))
        if not validation['passed']:
            self.assertFalse((out/'amy_statistics.tsv').exists())
            self.assertIn('STOPPED', json.loads((out/'amy_status.json').read_text())['status'])
        matrix = pd.read_csv(out/'confusion.tsv', sep='\t', index_col=0)
        self.assertEqual(int(matrix.to_numpy().sum()), validation['heldout_evaluable'])


if __name__ == '__main__':
    unittest.main()
