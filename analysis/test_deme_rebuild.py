"""Source-named regressions for encoding, path constraints and process covariance."""
import unittest

import numpy as np
import pandas as pd

from analysis.load_aadr import GenotypeMatrix, BUNDLE, MAIN_PREFIX, HO_PREFIX, parse_ind, iter_snp
from analysis.deme_rebuild import build_paths
from analysis.directional_consistency import slope_weights, process_variance, joint_score, conjunction_score


class TgenoEncodingTests(unittest.TestCase):
    def test_high_bits_first_and_arbitrary_site_order(self):
        packed = np.array([0b00011011,0b11100100],dtype=np.uint8)
        expected = np.array([0,1,2,3,3,2,1,0],dtype=np.uint8)
        np.testing.assert_array_equal(GenotypeMatrix._decode_bytes(packed),expected)
        sites = np.array([6,1,7,0,3])
        np.testing.assert_array_equal(GenotypeMatrix._decode_bytes(packed,sites),expected[sites])

    def test_real_cross_bundle_identical_calls(self):
        # The same individuals and SNP IDs occur in independently packed site orders.
        gm = GenotypeMatrix(f'{BUNDLE}/{MAIN_PREFIX}.geno')
        ho = GenotypeMatrix(f'{BUNDLE}/{HO_PREFIX}.geno')
        hi = {r.id:i for i,r in enumerate(parse_ind(f'{BUNDLE}/{HO_PREFIX}.ind'))}
        hs = {r[0]:i for i,r in enumerate(iter_snp(f'{BUNDLE}/{HO_PREFIX}.snp'))}
        left, right = [], []
        for i,r in enumerate(iter_snp(f'{BUNDLE}/{MAIN_PREFIX}.snp')):
            if r[0] in hs: left.append(i);right.append(hs[r[0]])
            if len(left)==1000: break
        checked = compared = 0
        for i,r in enumerate(parse_ind(f'{BUNDLE}/{MAIN_PREFIX}.ind')):
            if r.id not in hi: continue
            a,b = gm.individual_genotypes(i,np.array(left)),ho.individual_genotypes(hi[r.id],np.array(right))
            present = (a!=3)&(b!=3)
            np.testing.assert_array_equal(a[present],b[present])
            compared += int(present.sum());checked += 1
            if checked==20:break
        self.assertEqual(checked,20)
        self.assertGreater(compared,5000)


class AnchoredPathTests(unittest.TestCase):
    @staticmethod
    def fixture(bins,cells=None,centers=None):
        cells = cells or [0]*len(bins)
        centers = centers or [0]*len(bins)
        rows, pcs = [], []
        for b,cell,center in zip(bins,cells,centers):
            for i in range(10):
                rows.append(dict(lat_cell=cell,lon_cell=0,time_bin=b,date_BP=1000*b+50*i+100))
                pcs.append([center+(-1 if i%2 else 1),0])
        return pd.DataFrame(rows),np.array(pcs,dtype=float)

    def test_adjacent_bins_form_real_three_node_paths(self):
        meta,pcs = self.fixture([1,2,3])
        _,paths = build_paths(meta,pcs)
        self.assertEqual([len(p) for p in paths],[3])
        members = np.concatenate([n['members'] for n in paths[0]])
        self.assertEqual(len(np.unique(members)),30)
        self.assertGreaterEqual(meta.loc[members,'date_BP'].max()-meta.loc[members,'date_BP'].min(),2000)

    def test_missing_bin_and_geography_never_merge(self):
        meta,pcs = self.fixture([1,3,3],[0,0,1])
        _,paths = build_paths(meta,pcs)
        self.assertEqual([len(p) for p in paths],[1,1,1])

    def test_anchor_prevents_transitive_genetic_chain(self):
        meta,pcs = self.fixture([1,2,3,4],centers=[0,1,2,3])
        _,paths = build_paths(meta,pcs)
        self.assertEqual(sorted(len(p) for p in paths),[1,3])

    def test_duration_cap_prevents_giant_blob(self):
        meta,pcs = self.fixture(list(range(12)))
        _,paths = build_paths(meta,pcs)
        for path in paths:
            dates = meta.loc[np.concatenate([n['members'] for n in path]),'date_BP']
            self.assertLessEqual(dates.max()-dates.min(),8000)
        self.assertGreater(len(paths),1)


class DriftProcessTests(unittest.TestCase):
    def test_full_random_walk_covariance_not_independent_bins(self):
        time = np.array([0.,1.,3.,4.])
        weights = slope_weights(time,np.array([10,20,10,30]))
        sampling = np.array([.02,.01,.02,.01])
        q = np.array([.01,.04,.03])
        L = np.tril(np.ones((4,3)),k=-1)
        covariance = L @ np.diag(q) @ L.T + np.diag(sampling)
        self.assertAlmostEqual(process_variance(weights,sampling,q),weights @ covariance @ weights)
        self.assertAlmostEqual(weights.sum(),0)
        self.assertAlmostEqual(weights @ time,1)

    def test_forward_time_copy_number_direction(self):
        time = -np.array([6000,4000,2000])/1000
        weights = slope_weights(time,np.array([10,10,10]))
        self.assertGreater(weights @ np.array([3,5,7]),0)
        self.assertAlmostEqual(weights @ np.array([5,5,5]),0)

    def test_conjunction_cannot_be_driven_by_one_strong_deme(self):
        self.assertLess(conjunction_score([100.,-.1,-.2]),0)
        self.assertAlmostEqual(conjunction_score([100.,.1,.2]),.1)
        self.assertAlmostEqual(conjunction_score([-100.,-.1,-.2]),.1)

    def test_common_neutral_donor_can_force_parallel_direction(self):
        # Counterexample to 'arbitrary demography implies independent fair signs'.
        # Both resident sources and the donor are temporally constant.
        donor = .8
        mixture = np.array([0.,.2,.4])
        time = np.array([0.,1.,2.])
        weights = slope_weights(time,np.ones(3))
        for resident in [.1,.2,.3,.4,.5,.6]:
            path = (1-mixture)*resident + mixture*donor
            self.assertGreater(weights @ path,0)

    def test_shared_innovations_discount_joint_score(self):
        self.assertAlmostEqual(joint_score(np.ones(3),np.eye(3)),np.sqrt(3))
        self.assertAlmostEqual(joint_score(np.ones(3),np.ones((3,3))),1)


if __name__=='__main__':
    unittest.main()
