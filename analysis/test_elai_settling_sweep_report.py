"""Regression tests named after the trace-only settling analyser."""
import unittest

from elai_settling_sweep_report import CONFIGS, decision, parse_trace


def trace_from_deltas(deltas):
    likelihood = -1_000_000.0
    lines = [f' ### Exact update 0 0 {likelihood:.3f} {likelihood:.3f}']
    for i, delta in enumerate(deltas, 1):
        likelihood += delta
        lines.append(f' ### Exact update 0 {i} {likelihood:.3f} {delta:.3f}')
    return parse_trace('\n'.join(lines))


class ElaiSettlingSweepReportTests(unittest.TestCase):
    def test_parse_trace_skips_warm_and_first_bogus_increment(self):
        text = (' ### Approximate update0 0 -900 999\n'
                ' ### Exact update 0 0 -1e3 -1e3\n'
                ' ### Exact update 0 1 -9.99e2 1.0\n')
        trace = parse_trace(text)
        self.assertEqual([t['step'] for t in trace], [1, 2])
        self.assertIsNone(trace[0]['relative_increment'])
        self.assertEqual(trace[1]['relative_increment'], .001)
        self.assertEqual(trace[1]['absolute_increment'], 1)

    def test_decision_settles_decaying_nonnegative_small_final_window(self):
        stats = decision(trace_from_deltas([1000, 800, 300, 200, 50, 20, 10, 5]), True)
        self.assertTrue(stats['settles'])
        self.assertEqual(stats['quarter_size'], 2)
        self.assertEqual(stats['first_qualifying_step'], 8)
        self.assertEqual(stats['second_half_mean_abs_increment'], 21.25)
        self.assertEqual(stats['second_half_negative_fraction'], 0)

    def test_decision_constant_magnitude_plateau_does_not_decay(self):
        stats = decision(trace_from_deltas([10] * 12), True)
        self.assertTrue(stats['final_three_window_met'])
        self.assertFalse(stats['decays'])
        self.assertFalse(stats['settles'])

    def test_decision_oscillation_fails_even_with_final_window(self):
        stats = decision(trace_from_deltas([1000] * 8 + [-10, 10] * 3 + [5] * 3), True)
        self.assertTrue(stats['decays'])
        self.assertTrue(stats['final_three_window_met'])
        self.assertGreater(stats['second_half_negative_fraction'], .05)
        self.assertFalse(stats['settles'])

    def test_decision_negative_limit_is_strict(self):
        stats = decision(trace_from_deltas([1000] * 20 + [-1] + [1] * 19), True)
        self.assertEqual(stats['second_half_negative_fraction'], .05)
        self.assertTrue(stats['decays'])
        self.assertTrue(stats['final_three_window_met'])
        self.assertFalse(stats['settles'])

    def test_decision_early_hit_does_not_replace_final_window(self):
        stats = decision(trace_from_deltas([1000] * 4 + [5] * 7 + [200]), True)
        self.assertIsNotNone(stats['first_qualifying_step'])
        self.assertFalse(stats['final_three_window_met'])
        self.assertFalse(stats['settles'])

    def test_decision_incomplete_and_short_traces_cannot_settle(self):
        trace = trace_from_deltas([1000, 800, 300, 200, 50, 20, 10, 5])
        self.assertFalse(decision(trace, False)['settles'])
        for size in range(4):
            self.assertFalse(decision(trace[:size], True)['settles'])
        self.assertIsNone(decision([], False)['second_half_negative_fraction'])

    def test_parse_trace_rejects_gap_repetition_multirun_nan_malformed(self):
        bad = ['### Exact update 0 1 -100 1', '### Exact update 1 0 -100 1',
               '### Exact update 0 0 nan 0', '### Exact update 0 0 -100 nope',
               '### Exact update 0 0 -100 0\n### Exact update 0 0 -100 0']
        for text in bad:
            with self.subTest(text=text), self.assertRaises(ValueError):
                parse_trace(text)

    def test_plan_is_bounded_and_contains_requested_contrasts(self):
        self.assertEqual(CONFIGS[0], ('baseline', 1, 5, 0, 10, 70))
        self.assertGreaterEqual(sum(c[1] > 1 for c in CONFIGS), 2)
        self.assertIn((1, 3), [(c[1], c[2]) for c in CONFIGS])
        self.assertIn((1, 4), [(c[1], c[2]) for c in CONFIGS])
        self.assertTrue(any(c[3] > 0 for c in CONFIGS))
        self.assertTrue(any(c[4] != 10 for c in CONFIGS))
        for name, upper, lower, warm, _mg, steps in CONFIGS:
            with self.subTest(name=name):
                self.assertLessEqual((steps + warm) * 70 * ((upper * lower) / 5)**2 + 300, 5400)


if __name__ == '__main__':
    unittest.main()
