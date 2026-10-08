"""Regression tests for the registered one-hot residual implementation."""
import unittest
import numpy as np
from mixture_null import design, fit, simulation


class MeasuredMixtureTest(unittest.TestCase):
    def setUp(self):
        self.groups=np.repeat(['a','b','c'],50)
        self.t=np.tile(np.linspace(-4,0,50),3)+np.repeat([-3,0,3],50)
        self.g,self.n,self.x=design(self.t,self.groups)
        self.y=np.repeat([2,5,9],50)+.4*self.x+np.sin(np.arange(150))

    def test_measured_baseline_invariance(self):
        f=fit(self.t,self.y,self.groups)
        changed=fit(self.t,self.y+np.repeat([100,-200,300],50),self.groups)
        self.assertAlmostEqual(f['beta'],changed['beta'],places=12)
        self.assertAlmostEqual(f['se'],changed['se'],places=12)
        # Half-member means have exactly the same slope numerator.
        means=np.array([self.y[self.g==j][:25].mean() for j in range(3)])
        b=self.x@(self.y-means[self.g])/(self.x@self.x)
        self.assertAlmostEqual(b,f['beta'],places=12)

    def test_matches_explicit_fixed_effects(self):
        W=np.eye(3)[self.g]
        X=np.column_stack([W,self.t])
        coef=np.linalg.lstsq(X,self.y,rcond=None)[0]
        self.assertAlmostEqual(coef[-1],fit(self.t,self.y,self.groups)['beta'],places=12)
        f=fit(self.t,self.y,self.groups)
        self.assertEqual(f['df'],146)
        self.assertEqual(f['k'],3)

    def test_singleton_has_no_information(self):
        f=fit(self.t,self.y,self.groups)
        s=fit(np.r_[self.t,100],np.r_[self.y,10000],np.concatenate([self.groups,['singleton']]))
        self.assertAlmostEqual(f['beta'],s['beta'],places=12)
        self.assertAlmostEqual(f['se'],s['se'],places=12)
        self.assertEqual(f['df'],s['df'])
        self.assertFalse(fit(self.t,self.y,np.arange(150).astype(str))['identified'])

    def test_missing_values_are_not_alleles(self):
        y=self.y.copy(); y[0]=np.nan
        f=fit(self.t,y,self.groups)
        self.assertEqual(f['n'],149)
        self.assertTrue(f['identified'])

    def test_simulation_reproducible_and_falsifiable(self):
        # Same registered conditional simulation, toy design with equal sigma.
        a=simulation(self.t,self.groups,[2,5,9],[1,1,1],1,reps=100,seed=20261008)
        b=simulation(self.t,self.groups,[2,5,9],[1,1,1],1,reps=100,seed=20261008)
        self.assertEqual(a,b)
        self.assertGreater(a['power'],.8)
        z=simulation(self.t,self.groups,[2,5,9],[1,1,1],0,reps=1000,seed=20261008)
        self.assertLess(z['power'],.1)
        self.assertGreater(z['pooled_rejection'],.9)


if __name__=='__main__':
    unittest.main()
