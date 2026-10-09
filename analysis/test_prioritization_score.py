"""[observation] Software regressions for the registered prioritisation scorer."""
import unittest
import numpy as np
import pandas as pd
from mixture_null import fit
from prioritization_score import score_matrix, shape_path, rank_table


class PrioritizationScoreTests(unittest.TestCase):
    def setUp(self):
        self.meta=pd.DataFrame({'deme_id':np.repeat(np.arange(3),60),
            'node_id':np.tile(np.repeat(np.arange(3),20),3),
            't':np.tile(np.linspace(-5,-1,60),3)})
        self.frozen={'tau2':.00015,'directional_amplitude':.01,'R':np.eye(3)}
        rng=np.random.default_rng(24)
        self.Y=np.clip(.4+.035*(self.meta.t.to_numpy()[:,None]+3)+rng.normal(0,.06,(180,4)),.01,.99)

    def test_centered_hc3_matches_measured_mixture_fit(self):
        self.Y[0,0]=np.nan
        scored=score_matrix(self.Y,self.meta,self.frozen)
        for j in range(4):
            f=fit(self.meta.t,self.Y[:,j],self.meta.deme_id)
            self.assertAlmostEqual(scored.iloc[j].beta,f['beta'],places=12)
            self.assertAlmostEqual(scored.iloc[j].hc3_se,f['se'],places=12)

    def test_complement_invariant_ranking(self):
        a=score_matrix(self.Y,self.meta,self.frozen)
        b=score_matrix(1-self.Y,self.meta,self.frozen)
        for col in ['total_score','shape','consistency','N_eff','carriage_within_fraction']:
            np.testing.assert_allclose(a[col],b[col],atol=1e-10)
        np.testing.assert_allclose(a.beta,-b.beta,atol=1e-12)

    def test_smooth_ramp_and_step_are_distinguished(self):
        t=np.array([0.,1.,3.,4.])[:,None]
        ramp=t*.1
        step=np.array([0.,0.,.4,.4])[:,None]
        self.assertAlmostEqual(shape_path(t,ramp)[0][0],1.)
        self.assertAlmostEqual(shape_path(t,step)[0][0],0.)

    def test_partial_agreement_graded_not_unanimity(self):
        Y=self.Y.copy(); Y[120:]=1-Y[120:]
        a=score_matrix(Y,self.meta,self.frozen)
        self.assertTrue(a.valid.all())
        self.assertTrue((a.consistency>.55).all())
        self.assertTrue((a.consistency<.85).all())
        self.assertTrue(np.isfinite(a.total_score).all())

    def test_empty_background_path_does_not_create_observations(self):
        frozen={**self.frozen,'R':np.eye(4)}
        a=score_matrix(self.Y,self.meta,frozen)
        self.assertTrue(a.valid.all())
        self.assertTrue((a.n_lineages==3).all())
        self.assertTrue((a.n==180).all())
        np.testing.assert_allclose(a.N_eff,3.)
        np.testing.assert_allclose(a.breadth,.75)

    def test_missing_nodes_do_not_imply_a_jump_or_unanimity_filter(self):
        Y=self.Y.copy(); Y[self.meta.node_id.eq(1),0]=np.nan
        a=score_matrix(Y,self.meta,self.frozen)
        self.assertTrue(a.iloc[0].valid)
        self.assertEqual(a.iloc[0].n_shape_lineages,0)
        self.assertEqual(a.iloc[0]['shape'],0)
        self.assertTrue(np.isfinite(a.iloc[0].total_score))

    def test_missing_negative_support_and_carriage_never_excludes(self):
        Y=self.Y.copy(); Y[:120,0]=-9
        a=score_matrix(Y,self.meta,self.frozen)
        self.assertFalse(a.iloc[0].valid)
        self.assertEqual(a.iloc[0].total_score,-np.inf)
        self.assertTrue(a.valid.iloc[1:].all())
        self.assertTrue(np.isfinite(a.total_score.iloc[1:]).all())
        self.assertTrue((a.score_no_carriage.iloc[1:]>=a.score_carriage_half.iloc[1:]).all())
        self.assertTrue((a.score_carriage_half.iloc[1:]>=a.total_score.iloc[1:]).all())

if __name__=='__main__': unittest.main()
