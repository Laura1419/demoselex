#!/usr/bin/env python3
"""Artifact acceptance checks; no git probes and no raw-data mutation."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from two_channel import ROOT,OUT,inclusive_score


def main():
    d=pd.read_csv(ROOT/'data/derived/two_channel.tsv',sep='\t')
    source=pd.concat([pd.read_csv(f,sep='\t',usecols=['locus']) for prefix in ['shard','moderate'] for f in (ROOT/'data/derived/prio_panels').glob(prefix+'*_out/priority_loci.tsv')])
    assert len(d)==538262
    assert not d.locus.duplicated().any()
    assert set(d.locus)==set(source.locus)
    for col in ['beta_within','beta_between','beta_total','ordering_statistic','ordering_p','inclusive_score']:
        assert np.isfinite(d[col]).all(),col
    np.testing.assert_allclose(d.beta_total,d.beta_within+d.beta_between,atol=1e-12)
    np.testing.assert_allclose(d.inclusive_score,inclusive_score(d.z_within,d.z_between,d.ordering_p),atol=1e-12)
    assert d.inclusive_score.is_monotonic_decreasing
    assert d.ordering_p.between(1/4097,1).all()
    assert d.ordering_statistic.between(-1-1e-12,1+1e-12).all()
    blocks=pd.read_csv(OUT/'blocks.tsv',sep='\t')
    assert not blocks.duplicated(['chrom','block']).any()
    assert len(blocks)==10634
    assert blocks.control_only_block.sum()==626
    control=d.panel_class.eq('control')
    assert control.sum()==149994
    assert d.panel_class.eq('strong').sum()==31939
    assert d.panel_class.eq('moderate').sum()==356329
    ordering=binomtest(int((d.ordering_p.le(.05)&control).sum()),int(control.sum())).proportion_ci()
    bcontrol=blocks.control_only_block
    border=binomtest(int((blocks.ordering_p.le(.05)&bcontrol).sum()),int(bcontrol.sum())).proportion_ci()
    amy=json.loads((OUT/'amy_summary.json').read_text())
    assert abs(amy['beta_total']-amy['beta_within']-amy['beta_between'])<1e-12
    assert abs(amy['beta_pooled']-amy['beta_between']-amy['beta_within_pooled_contribution'])<1e-12
    assert amy['ordering']['ordering_donors']==7
    audit=json.loads((OUT/'akbari_audit.json').read_text())
    assert audit['AMY1']['n']==2204 and audit['AMY1']['strong_n']==0
    assert abs(audit['AMY1']['max_abs_X']-3.795)<.0001
    ordered=d[d.ordering_p.le(.05)].copy()
    ordered['block']=ordered.pos//250000
    ordered.drop_duplicates(['chrom','block'])[['locus','chrom','block','beta_within','beta_between','beta_total','ordering_statistic','ordering_p','rank']].to_csv(OUT/'ordered_blocks.tsv',sep='\t',index=False)
    ordered_block_keys=set(zip(ordered.chrom,ordered.block))
    control_block_keys=set(zip(blocks.loc[bcontrol,'chrom'],blocks.loc[bcontrol,'block']))
    result={'status':'passed','unique_scored_loci':len(d),'blocks':len(blocks),'control_loci':int(control.sum()),'control_only_blocks':int(bcontrol.sum()),'control_ordering_95CI':[float(ordering.low),float(ordering.high)],'control_block_ordering_95CI':[float(border.low),float(border.high)],'any_ordered_blocks':int(d[d.ordering_p.le(.05)].assign(block=lambda p:p.pos//250000).drop_duplicates(['chrom','block']).shape[0]),'control_any_ordered_blocks':len(ordered_block_keys&control_block_keys),'top_decile_channel_classes':{k:{'count':int((d.top_decile&d.channel_class.eq(k)).sum()),'control_count':int((d.top_decile&d.channel_class.eq(k)&control).sum()),'control_pass_fraction':float((d.top_decile&d.channel_class.eq(k))[control].mean())} for k in ['carriage','within','both']},'checks':['exact merged source coverage','finite complete columns','non-cancelling inclusive score','additive total identity','descending deterministic ranking','250kb unique block representatives','AMY pooled scaling identity','independent released AMY statistics']}
    (OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
