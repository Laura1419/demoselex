"""Source-named tests for preregistered CN-as-SNP validation and null logic."""
import unittest
import numpy as np
import pandas as pd
from cn_as_snp import (empirical_p, match_samples, ordering_gate, fit_binomial,
                       fit_continuous, completion_status)


class CnAsSnpTests(unittest.TestCase):
    def test_ordering_gate_requires_both_canonical_controls_and_convergence(self):
        rows = [dict(locus=n, model='GRM_on', z_score=z, converged=True)
                for n, z in [('LCT', -4), ('SLC45A2', 3), ('AMY_real_SNP', 1)]]
        self.assertTrue(ordering_gate(rows))
        self.assertFalse(ordering_gate(rows[:-1]))
        rows[1]['converged'] = False
        self.assertFalse(ordering_gate(rows))

    def test_ordering_gate_requires_prespecified_effect_separation(self):
        rows = [dict(locus=n, model='GRM_on', z_score=z, converged=True)
                for n, z in [('LCT', 2), ('SLC45A2', 2), ('AMY_real_SNP', 1.7)]]
        self.assertFalse(ordering_gate(rows))
        rows[-1]['z_score'] = 1.6
        self.assertTrue(ordering_gate(rows))
        rows[0]['z_score'] = 1.99
        self.assertFalse(ordering_gate(rows))

    def test_completion_distinguishes_failed_fits_from_budget_stop(self):
        rows = [dict(converged=True, permutations_completed=165, permutations_failed=34),
                dict(converged=False)]
        self.assertEqual(completion_status(rows, 199), 'complete_with_fit_failures')
        rows[0]['permutations_failed'] = 33
        self.assertEqual(completion_status(rows, 199), 'partial_permutation_budget')
        self.assertEqual(completion_status([dict(converged=True, permutations_completed=199)], 199), 'complete')

    def test_empirical_p_two_sided_plus_one_and_missing_null(self):
        self.assertEqual(empirical_p(-2, [1, -2, 3, np.nan]), .75)
        self.assertIsNone(empirical_p(2, []))

    def test_match_unique_prefix_exact_preference_and_ambiguity(self):
        frame = pd.DataFrame({'sample': ['A', 'B', 'C', 'D']})
        ids = np.array(['A.SG', 'B.AG', 'B.SG', 'C', 'C.SG'], dtype=object)
        got = match_samples(frame, ids)
        self.assertEqual(got.match_method.tolist(), ['unique_prefix', 'ambiguous', 'exact', 'absent'])
        self.assertEqual(got.iloc[2].pgen_id, 'C')

    def test_binomial_missing_sentinel_not_an_observed_allele(self):
        t = np.tile(np.arange(6), 8).astype(float)
        y = np.tile([0, 1, 0, 2, 1, 2], 8).astype(float)
        y[0] = -9
        got = fit_binomial(y, t, np.zeros((48, 48)))
        self.assertEqual(got['n'], 47)
        self.assertTrue(np.isfinite(got['z_score']))

    def test_continuous_independent_fit_matches_ols_wald(self):
        rng = np.random.default_rng(8)
        t = np.arange(40, dtype=float)
        y = .05 * t + rng.normal(size=40)
        got = fit_continuous(y, t, np.zeros((40, 40)))
        x = np.column_stack([np.ones(40), t])
        beta = np.linalg.lstsq(x, y, rcond=None)[0]
        sigma2 = np.sum((y - x @ beta) ** 2) / 38
        se = np.sqrt(sigma2 * np.linalg.inv(x.T @ x)[1, 1])
        self.assertAlmostEqual(got['z_score'], beta[1] / se, places=4)


if __name__ == '__main__':
    unittest.main()
