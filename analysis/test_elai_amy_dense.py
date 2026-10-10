"""Source-named regressions for dense QC order and convergence fail-closed rules."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
from elai_amy_dense import convergence, qc_mask, retain_spaced, consume


class DenseInstrumentTests(unittest.TestCase):
    def test_qc_before_spacing_keeps_common_neighbour(self):
        rows = [(0,1,100,'rare','A','G'), (1,1,101,'common','A','G'), (2,1,150,'too_close','A','G'), (3,1,201,'next','A','G')]
        h = np.array([[0,0,0,0], [0,1,1,2], [0,1,1,2], [0,1,1,2]], dtype=np.int8)
        mask, _, _ = qc_mask(h)
        self.assertEqual(retain_spaced(rows,h,mask,{},100),[1,3])

    def test_all_missing_fails_qc(self):
        mask, maf, call = qc_mask(np.array([[-9,-9],[0,2]],dtype=np.int8))
        self.assertEqual(mask.tolist(),[False,True])
        self.assertEqual(call.tolist(),[0,1])
        self.assertEqual(maf.tolist(),[0,.5])

    def test_spacing_persists_across_batches(self):
        last={}
        r=[(0,1,100,'a','A','G')]
        self.assertEqual(retain_spaced(r,None,np.array([True]),last,100),[0])
        self.assertEqual(retain_spaced([(1,1,199,'b','A','G')],None,np.array([True]),last,100),[])
        self.assertEqual(retain_spaced([(2,1,200,'c','A','G')],None,np.array([True]),last,100),[0])

    def test_convergence_requires_complete_fit_and_three_terminal_increments(self):
        trace=[dict(step=i+1,likelihood=v) for i,v in enumerate([-100000.,-99999.,-99998.,-99997.])]
        self.assertTrue(convergence(trace,True)['pass_gate'])
        self.assertFalse(convergence(trace,False)['pass_gate'])
        self.assertFalse(convergence(trace[:3],True)['pass_gate'])

    def test_negative_and_later_large_increment_do_not_converge(self):
        for ll in ([-10000.,-10000.1,-10000.2,-10000.3],[-100000.,-99999.,-99998.,-99997.,-99000.]):
            self.assertFalse(convergence([dict(likelihood=v) for v in ll],True)['pass_gate'])

    def test_full_export_stream_completeness_and_alternate_marker_compositions(self):
        markers=[(0,1,100,'m0'),(1,1,200,'m1')]
        with tempfile.TemporaryDirectory(prefix='dense-test-') as d:
            path=Path(d)/'tensor.tsv'
            with path.open('w') as f:
                f.write('individual_index\tindividual_id\tlocus\tchrom\tposition\tk1\tk2\tp_ordered\tq_unordered\n')
                for m in range(2):
                    for a in range(2):
                        for b in range(2):
                            f.write(f'0\t\tm{m}\t1\t{100+100*m}\t{a}\t{b}\t0.25\t{.25 if a==b else .5}\n')
            result={}; consume(path,1,2,markers,result)
            self.assertTrue(result['complete'])
            self.assertEqual(result['rows'],8)
            self.assertEqual(result['composition'],[[.5,.5]])
            self.assertEqual(result['thin_composition'],[[.5,.5]])
            path.write_text('individual_index\t\n')
            result={}; consume(path,1,2,markers,result)
            self.assertFalse(result['complete'])


if __name__=='__main__':
    unittest.main()
