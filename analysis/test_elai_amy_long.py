"""Driver regressions and live small-fixture marginal/joint equivalence."""
import csv
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import elai_amy_dense as dense
import elai_amy_long as driver


class LongDriverTests(unittest.TestCase):
    def test_fit_default_command_unchanged(self):
        class StopBeforeLaunch(Exception):
            pass
        with tempfile.TemporaryDirectory(dir=driver.BASE) as d:
            work = Path(d)
            (work/'output').mkdir()
            expected = [dense.BINARY, '-g', str(work/'panel.inp'), '-p', '1', '-pos', str(work/'panel.pos'),
                        '-C', '1', '-c', '5', '-w', '0', '-s', str(dense.STEPS), '-mg', '10', '-R', '20261010',
                        '-o', 'panel_k5_seed20261010', '--pkk-site', '1:104121000', '--pkk']
            with patch.object(dense, 'OUT', work), patch.object(dense.threading, 'Thread'), patch.object(dense.subprocess, 'Popen', side_effect=StopBeforeLaunch) as popen:
                with self.assertRaises(StopBeforeLaunch):
                    dense.fit(work, {'fitted_n': 12}, [], 5, 20261010, 900)
                self.assertEqual(popen.call_args.args[0], expected)

    def test_marginal_copy_counts_are_halved(self):
        with tempfile.TemporaryDirectory(dir=driver.BASE) as d:
            p = Path(d)/'counts.txt'
            p.write_text('1.854 0.146 0.200 1.800\n')
            result = {}
            dense.consume_marginal(p, 1, 2, [0, 1], result)
            self.assertTrue(result['complete'])
            np.testing.assert_allclose(result['composition'], [[.5135, .4865]])
            np.testing.assert_allclose(result['thin_composition'], [[.927, .073]])

    def test_marginal_truncation_is_failure(self):
        with tempfile.TemporaryDirectory(dir=driver.BASE) as d:
            p = Path(d)/'counts.txt'
            p.write_text('1.0 1.0\n')
            result = {}
            dense.consume_marginal(p, 1, 2, [0, 1], result)
            self.assertFalse(result['complete'])
            self.assertIn('incomplete', result['error'])

    def test_aggregate_not_tested_is_partial_not_pass(self):
        self.assertEqual(driver.aggregate([{'status': 'PASS'}, {'status': 'NOT TESTED'}]), 'PARTIAL')
        self.assertEqual(driver.aggregate([{'status': 'FAIL'}, {'status': 'NOT TESTED'}]), 'FAIL')

    def test_converged_gate_metric_path_and_marginal_coverage(self):
        # Simulate convergence ONLY inside isolated test artifacts. Real 3-step
        # runs retain their FAIL convergence and never authorize label use.
        source = driver.BASE/'work_fixture'
        with tempfile.TemporaryDirectory(dir=driver.BASE) as d:
            work = Path(d)
            (work/'output').mkdir()
            for name in ('panel.json', 'samples.json', 'markers.json', 'independent_pcs.npy'):
                (work/name).write_bytes((source/name).read_bytes())
            for k, seed in ((5, 20261010), (5, 20261011), (3, 20261010)):
                prefix = f'panel_k{k}_seed{seed}'
                r = json.loads((source/(prefix+'.run.json')).read_text())
                r['convergence']['pass_gate'] = True
                r['convergence']['status'] = 'PASS'
                r.update(export_mode='marginal', full_export_complete=False,
                         full_export_rows=0, marginal_export={'complete': True})
                driver.save(work/(prefix+'.run.json'), r)
                for suffix in ('.composition.tsv', '.thin_composition.tsv', '.pkk-site.tsv'):
                    (work/(prefix+suffix)).write_bytes((source/(prefix+suffix)).read_bytes())
                (work/'output'/(prefix+'.snpinfo.txt')).write_bytes((source/'output'/(prefix+'.snpinfo.txt')).read_bytes())
            with patch.object(dense, 'OUT', work):
                g = driver.gates(work)
            self.assertIsInstance(g['gates']['known_structure']['country_ami'], float)
            self.assertIn(g['gates']['stability']['seed']['status'], ('PASS', 'FAIL'))
            self.assertIn(g['gates']['stability']['cross_k']['status'], ('PASS', 'FAIL'))
            self.assertEqual(g['gates']['coverage']['components']['full_ordered_export']['status'], 'NOT TESTED')
            self.assertFalse(g['biological_use_authorized'])

    def test_live_fixture_joint_pk1_equivalence(self):
        source = driver.BASE/'work_fixture'
        self.assertTrue((source/'panel.inp').exists(), 'run fixture build first')
        panel, samples, markers, _ = driver.load(source)
        n, m, k = len(samples), len(markers), 5
        with tempfile.TemporaryDirectory(dir=driver.BASE) as d:
            work = Path(d)
            (work/'output').mkdir()
            commands = []
            for prefix, flags in [('joint', ['--pkk', '--pk1']), ('marginal', ['--pk1']), ('site_only', [])]:
                cmd = [dense.BINARY, '-g', str(source/'panel.inp'), '-p', '1', '-pos', str(source/'panel.pos'),
                       '-C', '1', '-c', str(k), '-w', '0', '-s', '3', '-mg', '10', '-R', '20261010',
                       '-o', prefix, '--pkk-site', '1:104121000', *flags]
                proc = subprocess.run(cmd, cwd=work, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
                self.assertEqual(proc.returncode, 1, proc.stdout)
                self.assertIn('-elai: finished', proc.stdout)
                commands.append(cmd)
            p = np.zeros((n, m, k, k)); q = p.copy()
            with (work/'output/joint.pkk.tsv').open() as f:
                rows = list(csv.DictReader(f, delimiter='\t'))
            self.assertEqual(len(rows), n*m*k*k)
            for rownum, r in enumerate(rows):
                i, a, b = int(r['individual_index']), int(r['k1']), int(r['k2'])
                j = (rownum//(k*k)) % m
                p[i,j,a,b] = float(r['p_ordered'])
                q[i,j,a,b] = float(r['q_unordered'])
            counts = np.array([np.fromstring(l, sep=' ').reshape(m,k)/2
                               for l in (work/'output/marginal.pk1.txt').read_text().splitlines()])
            ordered = p.sum(3)
            # Unique unordered pairs contribute half a copy to each off-diagonal state.
            unordered = .5*q.sum(3) + .5*np.diagonal(q, axis1=2, axis2=3)
            error = float(np.abs(ordered-counts).max())
            self.assertLessEqual(error, .00025+1e-12)  # 3-decimal copy counts / 2
            np.testing.assert_allclose(ordered, unordered, atol=1e-12, rtol=0)
            for suffix in ('.em.txt', '.admix.txt', '.ps21.txt', '.snpinfo.txt', '.pk1.txt'):
                self.assertEqual((work/'output'/('joint'+suffix)).read_bytes(), (work/'output'/('marginal'+suffix)).read_bytes(), suffix)
            for prefix in ('marginal', 'site_only'):
                self.assertEqual((work/'output/joint.pkk-site.tsv').read_bytes(), (work/'output'/(prefix+'.pkk-site.tsv')).read_bytes())
            joint_comp = ordered.mean(1)
            saved_joint = np.loadtxt(source/'panel_k5_seed20261010.composition.tsv')
            saved_marginal = np.loadtxt(driver.BASE/'work_fixture_marginal/panel_k5_seed20261010.composition.tsv')
            np.testing.assert_allclose(saved_joint, joint_comp, atol=1e-12, rtol=0)
            np.testing.assert_allclose(saved_marginal, counts.mean(1), atol=1e-12, rtol=0)
            evidence = dict(individuals=n, markers=m, k=k, ordered_rows=len(rows), commands=commands,
                marginal_max_error=error, allowed_error=.00025+1e-12,
                q_marginal_max_error=float(np.abs(ordered-unordered).max()),
                ordered_norm_error=float(np.abs(p.sum((2,3))-1).max()),
                site_only_byte_identical=True, legacy_outputs_byte_identical=True,
                joint_driver_composition_max_error=float(np.abs(saved_joint-joint_comp).max()),
                marginal_driver_composition_max_error=float(np.abs(saved_marginal-counts.mean(1)).max()))
            driver.save(driver.BASE/'long_equivalence.json', evidence)


if __name__ == '__main__':
    unittest.main()
