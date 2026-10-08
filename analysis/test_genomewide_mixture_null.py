"""Source-named regressions: per-locus support and frozen TRAIN-only variance."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import genomewide_mixture_null as gm


class GenomewideMixtureTest(unittest.TestCase):
    def test_public_aadr_reader_returns_alt(self):
        from load_aadr import GenotypeMatrix, BUNDLE, MAIN_PREFIX
        matrix=GenotypeMatrix(f'{BUNDLE}/{MAIN_PREFIX}.geno')
        sites=np.array([0,1,2,3,19,57,600,5000,100000])
        for person in [0,1,5]:
            raw=(matrix.packed_row(person)[sites//4] >> (6-2*(sites%4)).astype(np.uint8)) & 3
            expected=np.where(raw<=2,2-raw,3).astype(np.uint8)
            np.testing.assert_array_equal(matrix.individual_genotypes(person,sites),expected)

    def test_support_and_alt_orientation(self):
        t=np.tile(np.linspace(-4,0,60),2)
        group=np.repeat(['a','b'],60)
        meta=pd.DataFrame(dict(t=t,group=group))
        x=t+2
        signal=np.repeat([.3,.7],60)+.05*x
        Y=np.column_stack([signal,np.full(120,np.nan),np.repeat(-4.5,120)])
        # Six observed people do not make a supported component.
        Y[:6,1]=.5
        sites=pd.DataFrame(dict(id=['ramp','unsupported','missing_sentinel'],chrom=['2','3','4']))
        result=gm.score_matrix(Y,meta,sites,'test')
        self.assertEqual(list(result.id),['ramp'])
        self.assertAlmostEqual(result.beta.iloc[0],.05,places=12)
        inverse=gm.score_matrix(1-Y,meta,sites,'test')
        self.assertAlmostEqual(inverse.beta.iloc[0],-.05,places=12)
        self.assertEqual(result.n.iloc[0],120)
        self.assertEqual(result.k.iloc[0],2)

    def fixture(self,path,heldout_beta=0):
        N=200
        df=pd.DataFrame(dict(beta=np.r_[np.repeat(.02,100),np.repeat(heldout_beta,100)],
                             se=np.repeat(.005,N),af=np.repeat(.5,N),
                             partition=np.repeat(['train','heldout'],100)))
        for label in ['demes_aadr1240k','demes_akbari_pgen']:
            df.to_csv(path/f'genomewide_{label}.tsv',sep='\t',index=False)

    def test_no_heldout_parameter_fit(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            p=Path(d)
            with patch.object(gm,'OUT',p):
                self.fixture(p,0)
                first=gm.calibrated_deme_controls()
                self.fixture(p,100)
                second=gm.calibrated_deme_controls()
            for a,b in zip(first,second):
                self.assertEqual(a['tau2'],b['tau2'])
                self.assertEqual(a['training_loci'],100)
                self.assertEqual(a['shared_demographic_variance_parameters'],1)
                self.assertEqual(a['tau2_bootstrap_quantiles'],b['tau2_bootstrap_quantiles'])
                self.assertNotEqual(a['heldout_rejects'],b['heldout_rejects'])
                self.assertGreater(a['tau2'],0)

    def test_zero_excess_variance_boundary(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parent) as d:
            p=Path(d); self.fixture(p)
            for label in ['demes_aadr1240k','demes_akbari_pgen']:
                file=p/f'genomewide_{label}.tsv'
                df=pd.read_csv(file,sep='\t'); df['beta']=0.; df.to_csv(file,sep='\t',index=False)
            with patch.object(gm,'OUT',p):
                result=gm.calibrated_deme_controls()
            self.assertTrue(all(r['tau2']==0 for r in result))


if __name__=='__main__': unittest.main()
