"""Source-named regressions for the ELAI AMY input/export gate helpers."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from elai_amy_build import CORE, SITE, excluded, write_panel
import elai_amy_build
import elai_amy_gate
from elai_amy_gate import consume_full, joint_errors, read_joint, same_k, unlinked_parameters


class ElaiAmyTests(unittest.TestCase):
    def test_core_is_inclusive_and_virtual_inside(self):
        self.assertTrue(excluded(CORE[0]))
        self.assertTrue(excluded(CORE[1]))
        self.assertTrue(excluded(SITE))
        self.assertFalse(excluded(CORE[0] - 1))
        self.assertFalse(excluded(CORE[1] + 1))

    def test_bimbam_hardcalls_alt_coding_and_diploid_header(self):
        with tempfile.TemporaryDirectory(prefix='elai-amy-test-') as d:
            original = elai_amy_build.OUT
            try:
                elai_amy_build.OUT = Path(d)
                write_panel('test', [(0, 1, 95000000, 'snp', 'G', 'A')],
                            np.array([[0, 1, 2, -9]], np.int8), ['a', 'b', 'c', 'd'])
                lines = (Path(d) / 'test.inp').read_text().splitlines()
                self.assertEqual(lines, ['4', '1', 'IND a b c d', 'snp GG GA AA ??'])
                self.assertEqual((Path(d) / 'test.pos').read_text(), 'snp 95000000 1\n')
            finally:
                elai_amy_build.OUT = original

    def test_copy_count_vs_probability_normalization(self):
        p = np.array([[[.2, .3], [.3, .2]]])
        q = np.array([[[.2, .6], [.6, .2]]])
        e = joint_errors(p, q)
        self.assertEqual(e['ordered_sum_error'], 0)
        self.assertEqual(e['q_symmetry_error'], 0)
        self.assertEqual(e['q_triangle_sum_error'], 0)
        self.assertAlmostEqual(q.sum(), 1.6)  # square q is NOT normalized.
        self.assertAlmostEqual((2*p.sum(2)).sum(), 2)  # pk1 counts.

    def test_same_k_matching_is_label_invariant(self):
        a = np.array([[.9, .1], [.8, .2], [.1, .9], [.2, .8]])
        p = np.array([np.outer(v, v) for v in a])
        metrics = same_k(a, a[:, ::-1], p, p[:, ::-1, ::-1])
        self.assertEqual(metrics, {'ari': 1.0, 'mean_l1': 0.0, 'virtual_argmax_agreement': 1.0})

    def test_unlinked_control_preserves_parameter_line_structure(self):
        with tempfile.TemporaryDirectory(prefix='elai-amy-test-') as d:
            source, dest = Path(d) / 'linked.em', Path(d) / 'unlinked.em'
            source.write_text('0.1 0.2\n0.3 0.4\n0.5 0.6\n\n0.7 0.8\n')
            original = unlinked_parameters(source, dest, 1)
            self.assertEqual(original, '0.6')
            self.assertEqual(dest.read_text(), '0.1 0.2\n0.3 0.4\n0.5 1.0\n\n0.7 0.8\n')

    def test_full_stream_checks_sample_marker_order_and_counts(self):
        with tempfile.TemporaryDirectory(prefix='elai-amy-test-') as d:
            original = elai_amy_gate.OUT
            try:
                elai_amy_gate.OUT = Path(d)
                (Path(d) / 'samples.tsv').write_text('sample\na\nb\n')
                (Path(d) / 'markers.tsv').write_text('id\nm0\nm1\n')
                path = Path(d) / 'stream.tsv'
                lines = ['individual_index\tindividual_id\tlocus\tchrom\tposition\tk1\tk2\tp_ordered\tq_unordered\n']
                for i in range(2):
                    for marker in ('m0', 'm1'):
                        for a in range(2):
                            for b in range(2):
                                lines.append(f'{i}\t\t{marker}\t1\t95000000\t{a}\t{b}\t0.25\t0.25\n')
                path.write_text(''.join(lines))
                result = {}
                consume_full(path, 2, 2, 2, result)
                self.assertTrue(result['complete'])
                self.assertEqual(result['rows'], 16)
                self.assertEqual(result['individual_rows'], [8, 8])
                self.assertEqual(result['norm_error'], 0)
                path.write_text(''.join(lines).replace('m0', 'wrong', 1))
                result = {}
                consume_full(path, 2, 2, 2, result)
                self.assertFalse(result['complete'])
            finally:
                elai_amy_gate.OUT = original

    def test_joint_reader_rejects_incomplete_output(self):
        with tempfile.TemporaryDirectory(prefix='elai-amy-test-') as d:
            path = Path(d) / 'site.tsv'
            path.write_text('individual_index\tindividual_id\tlocus\tchrom\tposition\tk1\tk2\tp_ordered\tq_unordered\n'
                            '0\ta\tvirtual:1:104121000\t1\t104121000\t0\t0\t1\t1\n')
            with self.assertRaisesRegex(AssertionError, 'incomplete'):
                read_joint(path, 1, 2)

    def test_joint_reader_rejects_duplicate_rows(self):
        with tempfile.TemporaryDirectory(prefix='elai-amy-test-') as d:
            path = Path(d) / 'site.tsv'
            row = '0\ta\tvirtual:1:104121000\t1\t104121000\t0\t0\t1\t1\n'
            path.write_text('individual_index\tindividual_id\tlocus\tchrom\tposition\tk1\tk2\tp_ordered\tq_unordered\n' + row*2)
            with self.assertRaisesRegex(AssertionError, 'duplicate'):
                read_joint(path, 1, 2)


if __name__ == '__main__':
    unittest.main()
