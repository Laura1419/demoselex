#!/usr/bin/env python3
"""Two inclusive channels; preregistered in docs/two-channel-prereg.md."""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'data/derived/two_channel'
SEED = 20261009
REPS = 4096


def psd_correlation(R):
    R = (np.asarray(R)+np.asarray(R).T)/2
    v, U = np.linalg.eigh(R)
    R = (U*np.maximum(v,1e-8))@U.T
    d = np.sqrt(np.diag(R))
    return R/d[:,None]/d[None,:]


def ancestry_proxy(meta, group_col, prefix):
    """Aggregate independently measured old-trajectory covariance, not locus means."""
    old = pd.read_csv(ROOT/'data/derived/trajectories.tsv',sep='\t')
    old['bio'] = old.individual.str.split('.').str[0]
    old = old[old.trajectory_id.ge(0)]
    # Ambiguous sequencing duplicates must not produce multiple ancestry votes.
    ambiguity = old.groupby('bio').trajectory_id.nunique()
    old = old[old.bio.isin(ambiguity[ambiguity.eq(1)].index)].drop_duplicates('bio')
    lookup = old.set_index('bio').trajectory_id.to_dict()
    base = pd.read_csv(ROOT/'data/derived/traj_relatedness.tsv',sep='\t',index_col=0)
    base.columns = base.columns.astype(int); base.index = base.index.astype(int)
    labels = sorted(meta[group_col].unique())
    W = np.zeros((len(labels),len(base)))
    coverage = []
    for j,label in enumerate(labels):
        part = meta[meta[group_col].eq(label)]
        bio = part.biological_id if 'biological_id' in part else part['sample']
        mapped = [lookup[b] for b in bio if b in lookup]
        for x in mapped: W[j,base.index.get_loc(x)] += 1
        if mapped: W[j] /= len(mapped)
        coverage.append(dict(component=str(label),n=len(part),mapped=len(mapped),fraction=len(mapped)/len(part)))
    R = W@base.to_numpy()@W.T
    missing = W.sum(axis=1)==0
    R[np.ix_(missing,missing)] = .75
    np.fill_diagonal(R,1.)
    R = psd_correlation(R)
    pd.DataFrame(R,index=labels,columns=labels).to_csv(OUT/f'{prefix}_relatedness.tsv',sep='\t')
    pd.DataFrame(coverage).to_csv(OUT/f'{prefix}_mapping.tsv',sep='\t',index=False)
    return labels,R,coverage


def pair_stat(mu,R):
    """Signed weighted Kendall concordance; columns are loci or null draws."""
    mu = np.asarray(mu,float)
    if mu.ndim==1: mu=mu[:,None]
    a,b = np.triu_indices(len(R),1)
    w = np.maximum(0,1-R[a,b])
    return w@np.sign(mu[b]-mu[a])/max(w.sum(),1e-12)


def ordering(mu,t,n,R,seed=SEED,reps=REPS):
    """Actual correlated zero-direction drift simulations, three registered scenarios."""
    order = np.argsort(t,kind='stable')
    mu=np.asarray(mu)[order]; t=np.asarray(t)[order]; n=np.asarray(n)[order]
    R=R[np.ix_(order,order)]
    C=pair_stat(mu,R)
    u=t-t.min()
    B=(np.minimum.outer(u,u)+1)/np.sqrt((u[:,None]+1)*(u[None,:]+1))
    ps=[]
    for lam in [0,.5,1]:
        cov=psd_correlation((1-lam)*R+lam*B+np.diag(1/n))
        rng=np.random.default_rng(seed)
        sim=np.linalg.cholesky(cov+np.eye(len(n))*1e-10)@rng.normal(size=(len(n),reps))
        null=np.sort(abs(pair_stat(sim,R)))
        p=(1+reps-np.searchsorted(null,abs(C)-1e-12,side='left'))/(reps+1)
        ps.append(p)
    ranks=stats.rankdata(mu,axis=0)
    centered=ranks-ranks.mean(axis=0)
    tr=stats.rankdata(t); tr-=tr.mean()
    rho=(tr@centered)/np.maximum(np.sqrt((tr@tr)*(centered*centered).sum(axis=0)),1e-15)
    steps=(np.diff(mu,axis=0)>0).mean(axis=0)
    return dict(ordering_statistic=C,ordering_rho=rho,ordering_step_fraction=steps,
                ordering_p=np.max(ps,axis=0),ordering_p_genetic=ps[0],
                ordering_p_half_shared=ps[1],ordering_p_shared_drift=ps[2],
                ordering_neff=np.full(len(C),len(n)**2/R.sum()),ordering_donors=np.full(len(C),len(n)))


def inclusive_score(zw,zb,p):
    return np.log1p(np.asarray(zw)**2)+np.log1p(np.asarray(zb)**2)+np.maximum(0,-np.log10(p)-np.log10(20))


def compute_means(Y,meta):
    """Same per-locus called-member support as prioritization_score.score_matrix."""
    labels=sorted(meta.deme_id.unique()); M=Y.shape[1]
    mu=[]; count=[]; times=[]; supported=[]
    for g in labels:
        ix=np.flatnonzero(meta.deme_id.to_numpy()==g)
        A=np.asarray(Y[ix],float); ok=np.isfinite(A)&(A>=0)&(A<=1)
        A=np.where(ok,A,0); T=meta.t.to_numpy()[ix,None]
        n=ok.sum(axis=0); den=np.maximum(n,1)
        mt=(ok*T).sum(axis=0)/den
        lo=np.where(ok,T,np.inf).min(axis=0); hi=np.where(ok,T,-np.inf).max(axis=0)
        dates=np.zeros(M,int)
        for date in np.unique(T): dates+=ok[T[:,0]==date].any(axis=0)
        supported.append((n>=8)&(dates>=3)&((hi-lo)>=2))
        mu.append(A.sum(axis=0)/den); count.append(n); times.append(mt)
    return np.array(mu),np.array(times),np.array(count),np.array(supported)


def snp_run(limit=None):
    files=sorted((ROOT/'data/derived/prio_panels').glob('shard*_out/priority_loci.tsv'))+sorted((ROOT/'data/derived/prio_panels').glob('moderate*_out/priority_loci.tsv'))
    records=[]
    for f in files:
        table=pd.read_csv(f,sep='\t')
        panel=pd.read_csv(f.parent/'panel.tsv',sep='\t')
        meta=pd.read_csv(f.parent/'sample_metadata.tsv',sep='\t')
        Y=np.load(f.parent/'genotypes.npy',mmap_mode='r')
        if limit: table=table.head(limit)
        labels,R,coverage=ancestry_proxy(meta,'deme_id','snp')
        positions=panel.reset_index().set_index('locus')['index']
        ordering_rows=[]
        for start in range(0,len(table),2000):
            part=table.iloc[start:start+2000]
            indices=positions.loc[part.locus].to_numpy()
            mu,t,n,support=compute_means(Y[:,indices],meta)
            result={k:np.full(len(part),np.nan) for k in ['ordering_statistic','ordering_rho','ordering_step_fraction','ordering_p','ordering_p_genetic','ordering_p_half_shared','ordering_p_shared_drift','ordering_neff','ordering_donors']}
            patterns,inverse=np.unique(support.T,axis=0,return_inverse=True)
            for z,pattern in enumerate(patterns):
                donors=np.flatnonzero(pattern); cols=np.flatnonzero(inverse==z)
                if len(donors)<2: continue
                tt=t[np.ix_(donors,cols)]; nn=n[np.ix_(donors,cols)]
                # Fully imputed pgen panels share dates/counts. If not, do exact
                # per-column simulations rather than assume exchangeable designs.
                designs=np.vstack([tt,nn]).T
                unique,di=np.unique(designs,axis=0,return_inverse=True)
                for dz,design in enumerate(unique):
                    c=cols[di==dz]
                    fit=ordering(mu[np.ix_(donors,c)],design[:len(donors)],design[len(donors):],R[np.ix_(donors,donors)])
                    for k,v in fit.items(): result[k][c]=v
            ordering_rows.append(pd.DataFrame(result,index=part.index))
        ordered=pd.concat(ordering_rows).sort_index()
        d=table[['locus','chrom','pos','akbari_X','beta','composition_slope','corrected_se','hc3_se','n','n_lineages']].copy()
        d=d.rename(columns={'beta':'beta_within','composition_slope':'beta_between'})
        d['beta_total']=d.beta_within+d.beta_between
        d=pd.concat([d,ordered],axis=1)
        d.to_csv(OUT/f'{f.parent.name}.tsv',sep='\t',index=False)
        records.append(d)
        print('[observation] completed '+f.parent.name+' '+str(len(d)),flush=True)
    d=pd.concat(records,ignore_index=True).drop_duplicates('locus')
    control=d.akbari_X.abs().lt(3.16)
    scale=float(np.sqrt(np.mean(d.loc[control,'beta_between']**2)))
    d['z_within']=abs(d.beta_within)/d.corrected_se
    d['z_between']=abs(d.beta_between)/scale
    d['inclusive_score']=inclusive_score(d.z_within,d.z_between,d.ordering_p)
    d['panel_class']=np.where(control,'control',np.where(d.akbari_X.abs().gt(5.45),'strong','moderate'))
    w=d.z_within.ge(2); b=d.z_between.ge(2)
    d['channel_class']=np.select([b&~w,w&~b,b&w],['carriage','within','both'],default='neither')
    d['sign_disagreement']=d.beta_within*d.beta_between<0
    d=d.sort_values(['inclusive_score','locus'],ascending=[False,True]).reset_index(drop=True)
    d['rank']=np.arange(1,len(d)+1)
    d['top_decile']=d['rank'].le(int(np.ceil(len(d)/10)))
    d.to_csv(ROOT/'data/derived/two_channel.tsv',sep='\t',index=False)
    d['block']=d.pos//250000
    blocks=d.drop_duplicates(['chrom','block']).copy()
    onlycontrol=d.groupby(['chrom','block']).panel_class.apply(lambda s:s.eq('control').all())
    blocks['control_only_block']=[bool(onlycontrol.loc[(c,b)]) for c,b in zip(blocks.chrom,blocks.block)]
    blocks['block_rank']=np.arange(1,len(blocks)+1)
    blocks['top_decile']=blocks.block_rank.le(int(np.ceil(len(blocks)/10)))
    blocks.to_csv(OUT/'blocks.tsv',sep='\t',index=False)
    def metrics(frame,mask):
        criteria={'carriage':frame.channel_class.eq('carriage'),'within':frame.channel_class.eq('within'),'both':frame.channel_class.eq('both'),'neither':frame.channel_class.eq('neither'),'ordering':frame.ordering_p.le(.05),'ordered_carriage':frame.channel_class.eq('carriage')&frame.ordering_p.le(.05),'ordered_large_between':frame.z_between.ge(2)&frame.ordering_p.le(.05),'sign_disagreement':frame.sign_disagreement,'top_decile':frame.top_decile,'top_decile_ordered_carriage':frame.top_decile&frame.channel_class.eq('carriage')&frame.ordering_p.le(.05)}
        return {k:dict(count=int(v.sum()),fraction=float(v.mean()),control_count=int(v[mask].sum()),control_n=int(mask.sum()),control_pass_fraction=float(v[mask].mean()) if mask.sum() else None) for k,v in criteria.items()}
    summary=dict(n=len(d),panels=d.panel_class.value_counts().to_dict(),between_control_RMS=scale,variant_metrics=metrics(d,d.panel_class.eq('control')),block_n=len(blocks),block_metrics=metrics(blocks,blocks.control_only_block),relatedness_mapping=coverage)
    for cls in ['strong','moderate','control']:
        part=d[d.panel_class.eq(cls)]
        summary[cls]={'n':len(part),'median_abs_within':float(part.beta_within.abs().median()),'median_abs_between':float(part.beta_between.abs().median()),'median_total':float(part.beta_total.median()),'ordered_fraction':float(part.ordering_p.le(.05).mean()),'control_ordered_fraction':float(d.loc[d.panel_class.eq('control'),'ordering_p'].le(.05).mean())}
    (OUT/'genomewide_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def amy_run():
    from mixture_null import load_amy,fit,design
    ancient,cn,hap=load_amy()
    full=fit(ancient.t,ancient.AMY1,ancient.p2)
    groups=[]
    for label,part in ancient.groupby('p2'):
        groups.append(dict(component=label,mean_CN=float(part.AMY1.mean()),mean_t=float(part.t.mean()),**fit(part.t,part.AMY1,part.p2)))
    components=pd.DataFrame(groups).sort_values('mean_t')
    components.to_csv(OUT/'amy_components.tsv',sep='\t',index=False)
    t=ancient.t.to_numpy(); y=ancient.AMY1.to_numpy(float)
    g,n,x=design(t,ancient.p2); mu=np.bincount(g,weights=y)/n
    xp=t-t.mean(); between=float(xp@mu[g]/(xp@xp)); pooled=float(xp@y/(xp@xp))
    labels,R,coverage=ancestry_proxy(ancient,'p2','amy')
    selected=[j for j,k in enumerate(labels) if ancient.p2.eq(k).sum()>=8]
    cc=components.set_index('component').loc[[labels[j] for j in selected]]
    ordered=ordering(cc.mean_CN.to_numpy()[:,None],cc.mean_t.to_numpy(),cc.n.to_numpy(),R[np.ix_(selected,selected)])
    ordered={k:float(v[0]) for k,v in ordered.items()}
    rng=np.random.default_rng(SEED+1); bs=[]
    strata=[np.flatnonzero(g==j) for j in range(len(n))]
    for _ in range(1000):
        ix=np.concatenate([rng.choice(s,len(s),replace=True) for s in strata]); gg=g[ix]; tt=t[ix]; yy=y[ix]
        m=np.bincount(gg,weights=yy)/n; xx=tt-tt.mean(); bs.append(float(xx@m[gg]/(xx@xx)))
    se=float(np.std(bs,ddof=1))
    result=dict(n=len(ancient),beta_within=full['beta'],within_HC3_se=full['se'],within_p=full['p'],beta_between=between,beta_between_bootstrap_se=se,beta_total=full['beta']+between,beta_within_pooled_contribution=full['beta']*(x@x)/(xp@xp),beta_pooled=pooled,ordering=ordered,mapping=coverage,seed=SEED)
    result['inclusive_score']=float(inclusive_score(abs(full['beta'])/full['se'],abs(between)/se,ordered['ordering_p']))
    modern=cn[cn.type.eq('Modern')]
    result['modern_CN']={p:dict(n=int(modern.p1.eq(p).sum()),mean=float(modern.loc[modern.p1.eq(p),'AMY1'].mean())) for p in ['PEL','Pima','Surui']}
    slots=pd.concat([hap.rename(columns={k:'haplotype'}).drop(columns=[o]) for k,o in [('H1','H2'),('H2','H1')]],ignore_index=True)
    slots=slots[slots.type.eq('Modern')]
    result['H5']={name:dict(slots=int(mask.sum()),frequency=float(slots.loc[mask,'haplotype'].eq('H5').mean())) for name,mask in [('East_Asia',slots.p2.eq('EA')),('WEA',slots.p2.eq('WEA')),('Europe_5pop',slots.p1.isin(['CEU','FIN','GBR','IBS','TSI'])),('PEL',slots.p1.eq('PEL')),('Pima',slots.p1.eq('Pima'))]}
    (OUT/'amy_summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def akbari_audit():
    from load_akbari import iter_chunks
    windows={'AMY1':(1,103500000,104500000),'LCT':(2,136300000,136900000),'FADS':(11,61300000,61700000),'SLC45A2':(5,33800000,34200000),'SLC24A5':(15,48200000,48600000)}
    output={k:dict(n=0,max_abs_X=0.,strong_n=0) for k in windows}
    for d in iter_chunks(chunksize=500000,usecols=['CHROM','POS','ID','X']):
        for name,(c,lo,hi) in windows.items():
            p=d[d.CHROM.eq(c)&d.POS.between(lo,hi)]
            if len(p):
                r=output[name]; r['n']+=len(p); r['strong_n']+=int(p.X.abs().gt(5.45).sum())
                idx=p.X.abs().idxmax(); value=float(abs(p.loc[idx,'X']))
                if value>r['max_abs_X']: r['max_abs_X']=value; r['max_locus']=p.loc[idx,'ID']
    (OUT/'akbari_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    return output


def provenance(mode,outputs):
    sys.path.insert(0,str(ROOT)); import scripts.provenance as p
    p.get_git_info=lambda: {'disabled':'operator forbids every git command'}
    rec=p.RunRecorder('two-channel-'+mode,command=['.venv/bin/python','analysis/two_channel.py',*sys.argv[1:]])
    rec.set_seed(SEED)
    for f in ['docs/two-channel-prereg.md','analysis/two_channel.py','data/derived/traj_relatedness.tsv']:
        rec.record_input(f)
    rec.record_parameter('null_draws_per_scenario',REPS)
    rec.record_parameter('covariance_scenarios',[0,.5,1])
    if mode=='genomewide':
        for directory in sorted((ROOT/'data/derived/prio_panels').glob('*_out')):
            if not directory.name.startswith(('shard','moderate')): continue
            for name in ['priority_loci.tsv','panel.tsv','sample_metadata.tsv']:
                rec.record_input(str((directory/name).relative_to(ROOT)))
            large=directory/'genotypes.npy'
            rec.manifest['inputs'][str(large.relative_to(ROOT))]={'size_bytes':large.stat().st_size,'hash':'omitted: bounded immutable-input I/O'}
        outputs.append(ROOT/'data/derived/two_channel.tsv')
    elif mode=='amy':
        for name in ['copy_number_table_for_supplement.tsv','amylase_diversity_project-v1.1.zip']:
            rec.record_input('data/raw/amy/'+name)
    elif mode=='audit':
        large=ROOT/'data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz'
        rec.manifest['inputs'][str(large.relative_to(ROOT))]={'size_bytes':large.stat().st_size,'hash':'omitted: bounded immutable-input I/O'}
    rec.record_input('data/derived/trajectories.tsv')
    for f in outputs: rec.record_output(str(f.relative_to(ROOT)))
    rec.save()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['amy','genomewide','audit'],required=True)
    parser.add_argument('--limit',type=int)
    args=parser.parse_args(); OUT.mkdir(parents=True,exist_ok=True)
    if args.mode=='amy': result=amy_run()
    elif args.mode=='audit': result=akbari_audit()
    else: result=snp_run(args.limit)
    provenance(args.mode,list(OUT.glob('*.json'))+list(OUT.glob('*.tsv')))
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
