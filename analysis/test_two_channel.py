"""Source-named software checks for the two-channel protocol."""
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from mixture_null import fit,design
from two_channel import ordering,pair_stat,inclusive_score,compute_means,SEED,ROOT


class TwoChannelTests(unittest.TestCase):
    def test_centered_within_and_pooled_identity(self):
        t=np.r_[np.linspace(-8,-6,30),np.linspace(-4,0,40)]
        g=np.r_[np.repeat('a',30),np.repeat('b',40)]
        gi,n,x=design(t,g); mu=np.array([2.,9.])
        y=mu[gi]+.4*x
        f=fit(t,y,g)
        self.assertAlmostEqual(f['beta'],.4)
        f2=fit(t,y+np.array([10.,-10.])[gi],g)
        self.assertAlmostEqual(f2['beta'],.4)
        xp=t-t.mean(); measured=np.bincount(gi,weights=y)/n
        between=xp@measured[gi]/(xp@xp); pooled=xp@y/(xp@xp)
        self.assertAlmostEqual(pooled,between+f['beta']*(x@x)/(xp@xp))
        self.assertNotAlmostEqual(pooled,between+f['beta'])

    def test_inclusive_score_never_uses_cancellation(self):
        baseline=inclusive_score(0.,0.,.5)
        within=inclusive_score(4.,0.,.5)
        carriage=inclusive_score(0.,4.,.5)
        both=inclusive_score(4.,4.,.5)
        self.assertGreater(within,baseline)
        self.assertEqual(within,carriage)
        self.assertGreater(both,within)
        self.assertGreater(inclusive_score(4.,4.,.001),both)

    def test_ordering_null_determinism_and_orientation(self):
        R=np.eye(7); n=np.full(7,20); t=np.arange(7.)
        mu=np.arange(7.)[:,None]
        a=ordering(mu,t,n,R,reps=1024)
        b=ordering(mu,t,n,R,reps=1024)
        c=ordering(-mu,t,n,R,reps=1024)
        for k in a: np.testing.assert_array_equal(a[k],b[k])
        self.assertEqual(float(a['ordering_statistic'][0]),1.)
        self.assertEqual(float(c['ordering_statistic'][0]),-1.)
        np.testing.assert_array_equal(a['ordering_p'],c['ordering_p'])
        self.assertEqual(float(a['ordering_rho'][0]),1.)
        self.assertTrue(0<float(a['ordering_p'][0])<=1)

    def test_nested_donors_discount_effective_number_and_pairs(self):
        R=np.full((4,4),.9); np.fill_diagonal(R,1)
        a=ordering(np.arange(4.)[:,None],np.arange(4.),np.full(4,20),R,reps=128)
        self.assertLess(float(a['ordering_neff'][0]),2)
        # An identical pair gives no extra sign weight.
        R2=np.eye(3); R2[0,1]=R2[1,0]=1
        self.assertAlmostEqual(float(pair_stat([0,10,5],R2)[0]),0.)

    def test_cached_real_snp_composition_and_within(self):
        directory=ROOT/'data/derived/prio_panels/shard0_out'
        panel=pd.read_csv(directory/'panel.tsv',sep='\t')
        cached=pd.read_csv(directory/'priority_loci.tsv',sep='\t').head(12)
        meta=pd.read_csv(directory/'sample_metadata.tsv',sep='\t')
        index=panel.reset_index().set_index('locus')['index']
        Y=np.load(directory/'genotypes.npy',mmap_mode='r')[:,index.loc[cached.locus].to_numpy()]
        mu,t,n,support=compute_means(Y,meta)
        labels=sorted(meta.deme_id.unique())
        for j,row in enumerate(cached.itertuples()):
            groups=np.array(labels)[support[:,j]]
            ok=meta.deme_id.isin(groups).to_numpy()&np.isfinite(Y[:,j])&(Y[:,j]>=0)&(Y[:,j]<=1)
            tt=meta.t.to_numpy()[ok]; yy=Y[ok,j]; gg=meta.deme_id.to_numpy()[ok]
            f=fit(tt,yy,gg)
            self.assertAlmostEqual(f['beta'],row.beta,places=10)
            means={g:yy[gg==g].mean() for g in groups}
            xp=tt-tt.mean(); between=xp@np.array([means[g] for g in gg])/(xp@xp)
            self.assertAlmostEqual(between,row.composition_slope,places=10)

if __name__=='__main__': unittest.main()
