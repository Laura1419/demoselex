"""Synthetic, source-named checks; none validate biological state assignment."""
import unittest
import numpy as np
from analysis.local_ancestral_state_checks import (
    genotype_probabilities, pair_marginal_distribution,
    shared_direction_statistic, spectral_effective_count,
    ideal_frequency_floor, outside_amy_exclusion, assert_amy_features_excluded,
)


class LocalAncestralStateAlgebraTests(unittest.TestCase):
    def test_diploid_convolution_matches_enumeration(self):
        for p_a in [0., .1, .5, .9, 1.]:
            for p_b in [0., .1, .5, .9, 1.]:
                exact = np.zeros(3)
                for a in [0, 1]:
                    for b in [0, 1]:
                        exact[a+b] += (p_a if a else 1-p_a)*(p_b if b else 1-p_b)
                np.testing.assert_allclose(genotype_probabilities(p_a, p_b), exact)

    def test_mixed_only_named_frequency_not_identifiable(self):
        np.testing.assert_allclose(genotype_probabilities(.9, .1), [.09, .82, .09])
        np.testing.assert_allclose(genotype_probabilities(.9, .1),
                                   genotype_probabilities(.1, .9))
        # In an AB-only cohort the named A direction can be reversed without
        # changing either early or late observable genotype distributions.
        np.testing.assert_allclose(genotype_probabilities(.8, .2),
                                   genotype_probabilities(.2, .8))
        self.assertGreater(.9-.8, 0)
        self.assertLess(.1-.2, 0)

    def test_homo_state_anchors_break_swap_symmetry(self):
        self.assertFalse(np.allclose(genotype_probabilities(.9, .9),
                                     genotype_probabilities(.1, .1)))

    def test_pair_dosages_do_not_determine_joint_distribution(self):
        # Both q distributions imply E[n_A]=E[n_B]=1, but very different G.
        mixed = pair_marginal_distribution([.9, .1], {(0, 1): 1.})
        pure = pair_marginal_distribution([.9, .1], {(0, 0): .5, (1, 1): .5})
        np.testing.assert_allclose(mixed, [.09, .82, .09])
        np.testing.assert_allclose(pure, [.41, .18, .41])
        self.assertFalse(np.allclose(mixed, pure))

    def test_naive_ancestry_weighting_is_not_asafe(self):
        # One chromosome from each state, true p_A=.9,p_B=.1: naive G/2
        # allocation would estimate .5 for BOTH states in expectation.
        g_prob = genotype_probabilities(.9, .1)
        naive = float(np.arange(3) @ g_prob / 2)
        self.assertAlmostEqual(naive, .5)
        self.assertNotAlmostEqual(naive, .9)

    def test_marginalize_pair_before_log(self):
        q = {(0, 0): .25, (0, 1): .50, (1, 1): .25}
        result = pair_marginal_distribution([.8, .2], q)
        np.testing.assert_allclose(result, [.25, .5, .25])
        self.assertNotAlmostEqual(np.log(result[1]),
            sum(w*np.log(genotype_probabilities([.8,.2][a], [.8,.2][b])[1])
                for (a,b),w in q.items()))
        for invalid in [{(0, 1): .5}, {(1, 0): 1.}, {(0, 1): -1., (1,1): 2.}]:
            with self.assertRaises(ValueError):
                pair_marginal_distribution([.8, .2], invalid)

    def test_common_direction_not_raw_country_count(self):
        self.assertAlmostEqual(shared_direction_statistic([.3,.1,.2]), .1)
        self.assertAlmostEqual(shared_direction_statistic([-.3,-.1]), .1)
        self.assertLess(shared_direction_statistic([-.1,.3]), 0)
        self.assertEqual(shared_direction_statistic([0,.1]), 0)
        with self.assertRaises(ValueError):
            shared_direction_statistic([.1])

    def test_spectral_count_independence_and_duplicate_states(self):
        self.assertAlmostEqual(spectral_effective_count(np.eye(5)), 5.)
        self.assertAlmostEqual(spectral_effective_count(np.ones((5,5))), 1.)
        self.assertAlmostEqual(spectral_effective_count([[1,.5],[.5,1]]), 1.6)
        # Negative correlation is not clipped away; dependence still counts.
        self.assertAlmostEqual(spectral_effective_count([[1,-.5],[-.5,1]]), 1.6)
        with self.assertRaises(ValueError):
            spectral_effective_count([[1,2],[2,1]])
        with self.assertRaises(ValueError):
            spectral_effective_count([])

    def test_floor_units_and_support_loss(self):
        baseline = ideal_frequency_floor(112,113)
        self.assertAlmostEqual(baseline, .1320697, places=5)
        self.assertGreater(ideal_frequency_floor(11,11), baseline)
        self.assertAlmostEqual(ideal_frequency_floor(448,452), baseline/2)
        with self.assertRaises(ValueError):
            ideal_frequency_floor(0,113)

    def test_runtime_assertion_rejects_corrected_core_and_boundaries(self):
        for position in [103915000, 104255001, 104301314, 104327000]:
            for chromosome in [1, '1', 'chr1']:
                with self.assertRaisesRegex(AssertionError, 'AMY exclusion violated'):
                    assert_amy_features_excluded([chromosome], [position])
        assert_amy_features_excluded(['chr1', '1', 'chr2'],
                                     [103914999, 104327001, 104301314])

    def test_runtime_assertion_rejects_unknown_or_malformed_coordinates(self):
        for chromosome in ['unknown', 'chrX', '1.0']:
            with self.assertRaises(ValueError):
                assert_amy_features_excluded([chromosome], [104301314])
        for position in [np.nan, np.inf, -1, 104301314.5]:
            with self.assertRaises(ValueError):
                assert_amy_features_excluded([1], [position])
        with self.assertRaises(ValueError):
            assert_amy_features_excluded([1, 2], [103914999])

    def test_corrected_inclusive_exclusion_and_perturbation(self):
        positions = np.array([103914999,103915000,104255001,104301314,104327000,104327001])
        keep = outside_amy_exclusion(np.ones(6), positions)
        np.testing.assert_array_equal(keep, [True,False,False,False,False,True])
        dosage = np.arange(12).reshape(6,2).copy()
        perturbed = dosage.copy(); perturbed[~keep] = 999
        np.testing.assert_array_equal(dosage[keep], perturbed[keep])
        # Same coordinate on another chromosome is not AMY.
        self.assertTrue(outside_amy_exclusion([2], [104301314])[0])


if __name__ == '__main__':
    unittest.main()
