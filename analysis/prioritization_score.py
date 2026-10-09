#!/usr/bin/env python3
"""[inference] Frozen ranking protocol: docs/prioritization-prereg.md.

[observation] No significance filter or locus-fitted evolutionary parameters.
[inference] The LLR is a common-slope approximation, not a jump/path likelihood.
"""
from __future__ import annotations
import argparse
import json
import os
import pickle
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import ndtr
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/derived/prioritization'
SEED = 20261009
CONTROLS = [
    ('LCT_MCM6',2,136300000,136900000,136608646),
    ('SLC24A5',15,48200000,48600000,48426484),
    ('SLC45A2',5,33800000,34200000,33951693),
    ('TYR',11,88800000,89200000,89017961),
    ('OCA2_HERC2',15,28000000,28600000,28365618),
    ('FADS1_FADS2',11,61300000,61700000,61570783),
    ('TLR1_TLR6',4,38500000,39000000,38800000),
    ('AMY1_region',1,103500000,104500000,104000000)]


def frozen_inputs():
    """[observation] Read prior TRAIN law without refitting it on new candidates."""
    prior = json.loads((ROOT/'analysis/mixture_null/genomewide_corrected_results.json').read_text())
    law = next(r for r in prior if r['bundle']=='demes_akbari_pgen')
    train = pd.read_csv(ROOT/'analysis/mixture_null/genomewide_demes_akbari_pgen.tsv',sep='\t')
    train = train[train.partition.eq('train')]
    a = float(np.sqrt(np.mean(train.beta**2)))
    nodes = pd.read_csv(ROOT/'analysis/deme_rebuild/nodes.tsv',sep='\t')
    pc = nodes.groupby('deme_id')[[f'PC{i}' for i in range(1,9)]].mean().sort_index().to_numpy()
    dist2 = np.sum((pc[:,None,:]-pc[None,:,:])**2,axis=2)
    scale2 = np.median(dist2[dist2>0])
    R = np.exp(-dist2/(2*scale2))
    return dict(tau2=law['tau2'], directional_amplitude=a, training_loci=len(train), R=R)


def shape_path(times, values):
    """[inference] Time-proportional ramp=1, one step=0; no step fitting."""
    times, values = np.asarray(times), np.asarray(values)
    delta = np.diff(values,axis=0)
    total = np.sum(abs(delta),axis=0)
    monotonicity = np.divide(abs(values[-1]-values[0]),total,out=np.zeros_like(total),where=total>1e-15)
    concentration = np.divide(np.max(abs(delta),axis=0),total,out=np.ones_like(total),where=total>1e-15)
    ramp_share = np.max(np.diff(times,axis=0),axis=0)/(times[-1]-times[0])
    smoothness = np.clip((1-concentration)/np.maximum(1-ramp_share,1e-12),0,1)
    return monotonicity*smoothness, monotonicity, smoothness


def score_matrix(Y, meta, frozen):
    """[inference] Vectorised measured-baseline HC3, equivalent to mixture_null.fit.

    [observation] Missingness-specific means/time centring and support are recomputed.
    """
    Y = np.asarray(Y,float).copy()
    Y[(Y<0)|(Y>1)|~np.isfinite(Y)] = np.nan
    t = meta.t.to_numpy(float)
    group = meta.deme_id.to_numpy(int)
    K = frozen['R'].shape[0]
    M = Y.shape[1]
    numer = np.zeros(M); sxx = np.zeros(M); total_n = np.zeros(M)
    sy = np.zeros(M); st = np.zeros(M); stt = np.zeros(M); sty = np.zeros(M)
    stmu = np.zeros(M); groups = []; qs = []; supported = []
    shape_sum = np.zeros(M); mono_sum=np.zeros(M); smooth_sum=np.zeros(M); shape_n=np.zeros(M)
    time_lo=np.full(M,np.inf); time_hi=np.full(M,-np.inf)
    for g in range(K):
        ix = np.flatnonzero(group==g)
        if len(ix)==0:
            qs.append((np.zeros(M),np.zeros(M)))
            supported.append(np.zeros(M,bool))
            continue
        T = t[ix,None]; A=Y[ix]; ok=np.isfinite(A); A0=np.nan_to_num(A)
        n=ok.sum(axis=0); den=np.maximum(n,1)
        mt=(ok*T).sum(axis=0)/den; mu=A0.sum(axis=0)/den
        x=(T-mt)*ok
        xx=(x*x).sum(axis=0); xy=(x*A0).sum(axis=0)
        lo=np.min(np.where(ok,T,np.inf),axis=0); hi=np.max(np.where(ok,T,-np.inf),axis=0)
        dates=np.zeros(M,int)
        for date in np.unique(t[ix]): dates += ok[t[ix]==date].any(axis=0)
        support=(n>=8)&(dates>=3)&((hi-lo)>=2)&(xx>1e-12)
        beta=np.divide(xy,xx,out=np.zeros(M),where=xx>1e-12)
        e=(A0-mu-beta*x)*ok
        h=np.divide(1.,den)[None,:]+np.divide(x*x,np.maximum(xx,1e-12))
        adj=np.divide(x*e,1-h,out=np.zeros_like(e),where=(h<1-1e-12)&ok)
        se=np.sqrt((adj*adj).sum(axis=0))/np.maximum(xx,1e-12)
        qs.append((beta,se)); supported.append(support)
        numer+=xy*support; sxx+=xx*support; total_n+=n*support
        sy+=A0.sum(axis=0)*support; st+=(ok*T).sum(axis=0)*support
        stt+=(ok*T*T).sum(axis=0)*support; sty+=(A0*T).sum(axis=0)*support
        stmu+=((ok*T).sum(axis=0)*mu)*support
        time_lo=np.minimum(time_lo,np.where(support,lo,np.inf)); time_hi=np.maximum(time_hi,np.where(support,hi,-np.inf))
        groups.append((A0,ok,n,x,mu,support))
        # Nodes and node dates are measured using the same called members.
        vals=[]; times=[]; node_ok=[]
        for node in sorted(meta.iloc[ix].node_id.unique()):
            w=meta.iloc[ix].node_id.to_numpy()==node
            nn=ok[w].sum(axis=0); d=np.maximum(nn,1)
            vals.append(A0[w].sum(axis=0)/d)
            times.append((ok[w]*T[w]).sum(axis=0)/d)
            node_ok.append(nn>=3)
        if len(vals)>=3:
            vals=np.array(vals); times=np.array(times); node_ok=np.array(node_ok)
            # Pattern grouping retains all supported >=3-call nodes, not just complete paths.
            patterns, inv=np.unique(node_ok.T,axis=0,return_inverse=True)
            for p,pattern in enumerate(patterns):
                cols=np.flatnonzero(inv==p)
                rows=np.flatnonzero(pattern)
                if len(rows)<3: continue
                vv=vals[np.ix_(rows,cols)]; tt=times[np.ix_(rows,cols)]
                order=np.argsort(tt,axis=0)
                vv=np.take_along_axis(vv,order,axis=0); tt=np.take_along_axis(tt,order,axis=0)
                H,mo,sm=shape_path(tt,vv)
                valid=support[cols] & (np.ptp(tt,axis=0)>0)
                shape_sum[cols]+=H*valid; mono_sum[cols]+=mo*valid; smooth_sum[cols]+=sm*valid; shape_n[cols]+=valid
    support=np.array(supported); k=support.sum(axis=0)
    beta=numer/np.maximum(sxx,1e-12)
    hc3=np.zeros(M)
    for A0,ok,n,x,mu,s in groups:
        e=(A0-mu-beta*x)*ok
        h=1/np.maximum(n,1)[None,:]+x*x/np.maximum(sxx,1e-12)
        adj=np.divide(x*e,1-h,out=np.zeros_like(e),where=(h<1-1e-12)&ok)
        hc3+=(adj*adj).sum(axis=0)*s
    se=np.sqrt(hc3)/np.maximum(sxx,1e-12)
    af=sy/np.maximum(total_n,1)
    V=np.maximum(se*se+frozen['tau2']*af*(1-af),1e-16)
    a=frozen['directional_amplitude']
    llr=np.logaddexp(a*beta/V,-a*beta/V)-np.log(2)-a*a/(2*V)
    q=np.array([ndtr(np.sign(beta)*b/np.sqrt(np.maximum(s*s+frozen['tau2']*af*(1-af),1e-16))) for b,s in qs])
    consistency=(q*support).sum(axis=0)/np.maximum(k,1)
    denom=np.einsum('im,ij,jm->m',support.astype(float),frozen['R']**2,support.astype(float),optimize=True)
    neff=k*k/np.maximum(denom,1e-12); breadth=neff/K
    pooled_sxx=stt-st*st/np.maximum(total_n,1)
    comp=(stmu-st*sy/np.maximum(total_n,1))/np.maximum(pooled_sxx,1e-12)
    within=beta*sxx/np.maximum(pooled_sxx,1e-12)
    W=abs(within)/np.maximum(abs(within)+abs(comp),1e-15)
    H=shape_sum/np.maximum(shape_n,1)
    base=llr+np.log(.25+.75*consistency)+np.log(.5+.5*breadth)+np.log(.25+.75*H)
    valid=(total_n>=100)&(k>=2)&(af>0)&(af<1)&(sxx>0)
    score=base+np.log(.25+.75*W)
    score[~valid]=-np.inf
    return pd.DataFrame(dict(total_score=score, score_carriage_half=base+np.log(.5+.5*W),
        score_no_carriage=base, directional_llr=llr, magnitude_z=abs(beta)/np.sqrt(V),
        magnitude_slope=abs(beta), beta=beta, hc3_se=se, corrected_se=np.sqrt(V),
        consistency=consistency, breadth=breadth, shape=H, monotonicity=mono_sum/np.maximum(shape_n,1),
        smoothness=smooth_sum/np.maximum(shape_n,1), n_shape_lineages=shape_n.astype(int),
        carriage_within_fraction=W, composition_slope=comp, n_lineages=k, N_eff=neff,
        n=total_n.astype(int), time_span_kyr=time_hi-time_lo, af=af, valid=valid))


def verify_variant_map(panel):
    """[observation] Validate cached indices against immutable pvar allele IDs."""
    from load_akbari_geno import iter_pvar, paths
    expected=dict(zip(panel.pgen_index.astype(int),panel.locus))
    remaining=set(expected)
    for i,(chrom,pos,vid,ref,alt) in enumerate(iter_pvar(paths()[1])):
        if i in remaining:
            if vid!=expected[i] or vid!=f'{chrom}_{pos}_{ref}_{alt}':
                raise ValueError(f'cached pgen index/allele mismatch at {i}')
            remaining.remove(i)
            if not remaining: return len(expected)
    raise ValueError('cached pgen indices exceed pvar length')


def panel_and_data(panel_file=None):
    from load_akbari_geno import AkbariGeno
    from load_aadr import iter_snp, MAIN_PREFIX
    pool=pickle.load(open(ROOT/'data/derived/tier2_pool_cache.pkl','rb'))
    pool=pool[pool.pgen_index.ge(0)].copy()
    embedding=set(np.load(ROOT/'analysis/deme_rebuild/panels.npz')['embedding_sites'].tolist())
    exclude=set()
    for i,s in enumerate(iter_snp(str(ROOT/'data/raw/aadr-v66.p1'/f'{MAIN_PREFIX}.snp'))):
        if i in embedding: exclude.add(f'{s[1]}_{s[3]}_{s[4]}_{s[5]}')
    pool=pool[~pool.ID.isin(exclude)].copy()
    if panel_file is not None:
        # [inference] External caller supplies the variant panel (e.g. a shard of the
        # strong-significance set). Scoring math is unchanged; only the input set differs.
        panel=pd.read_csv(panel_file,sep='\t')
        need={'locus','chrom','pos','akbari_X','pgen_index'}
        if not need.issubset(panel.columns):
            raise ValueError(f'panel file must have {sorted(need)}')
        panel=panel[['locus','chrom','pos','akbari_X','pgen_index']].copy()
        panel['pgen_index']=panel.pgen_index.astype(int)
    else:
        pool['block']=pool.POS//250000
        pool['distance']=abs(pool.POS-(pool.block*250000+125000))
        chosen=pool.sort_values(['distance','POS','ID']).drop_duplicates(['CHROM','block'])
        anchors=[(c,a) for _,c,_,_,a in CONTROLS]+[(1,x) for x in [103625000,103875000,104125000,104375000]]
        extra=[]
        for c,a in anchors:
            d=pool[pool.CHROM.eq(c)]
            extra.append(d.loc[(d.POS-a).abs().idxmin()])
        panel=pd.concat([chosen,pd.DataFrame(extra)]).drop_duplicates('ID').sort_values('pgen_index').reset_index(drop=True)
        panel=panel.rename(columns={'ID':'locus','CHROM':'chrom','POS':'pos','X':'akbari_X'})
    OUT.mkdir(parents=True,exist_ok=True)
    panel.to_csv(OUT/'panel.tsv',sep='\t',index=False)
    meta=pd.read_csv(ROOT/'analysis/deme_rebuild/individuals.tsv',sep='\t')
    meta=meta[meta.deme_id.ge(0)].copy(); meta['deme_id']=meta.deme_id.astype(int); meta['t']=-meta.date_BP/1000
    if meta.individual.duplicated().any() or meta.biological_id.duplicated().any():
        raise ValueError('duplicate longitudinal individual identity')
    if not np.isfinite(meta.t).all() or not (meta.t<0).all():
        raise ValueError('invalid ancient date')
    geno=AkbariGeno(); index={s:i for i,s in enumerate(geno.sample_ids)}
    meta=meta[meta.individual.isin(index)].copy(); meta['pgen_row']=meta.individual.map(index)
    meta=meta.sort_values('pgen_row').reset_index(drop=True)
    geno.change_sample_subset(meta.pgen_row.to_numpy(np.uint32))
    Y=geno.read_variants(panel.pgen_index.to_numpy(np.uint32)).T.astype(float)/2
    geno.close(); Y[(Y<0)|(Y>1)]=np.nan
    meta.to_csv(OUT/'sample_metadata.tsv',sep='\t',index=False)
    np.save(OUT/'genotypes.npy',Y)
    return panel,meta,Y,dict(cache_variants=len(pool),cache_significant=int((pool.X.abs()>5.45).sum()),panel_variants=len(panel),shared_people=len(meta))


def ns_for(M): return sorted(set([n for n in [10,25,50,100,250,500,int(np.ceil(M/10))] if n<=M]))


def rank_table(panel,scores):
    table=pd.concat([panel[['locus','chrom','pos','akbari_X']],scores],axis=1)
    table=table[table.valid].copy().sort_values(['total_score','locus'],ascending=[False,True])
    table['rank']=np.arange(1,len(table)+1)
    table['percentile']=100*(1-(table['rank']-1)/len(table))
    table['akbari_abs_X']=table.akbari_X.abs()
    for col in ['score_carriage_half','score_no_carriage']:
        table['rank_'+col]=table[col].rank(method='first',ascending=False).astype(int)
    return table


def validate_controls(table):
    hits=np.zeros(len(table),bool); regions=[]; recall=[]
    for name,c,lo,hi,_ in CONTROLS:
        d=table[table.chrom.eq(c)&table.pos.between(lo,hi)]
        hits |= (table.chrom.eq(c)&table.pos.between(lo,hi)).to_numpy()
        regions.append(dict(control=name,variants=len(d),best_locus=None if d.empty else d.iloc[0].locus,
            rank=None if d.empty else int(d.iloc[0]['rank']),percentile=None if d.empty else float(d.iloc[0].percentile)))
    for n in ns_for(len(table)):
        r=sum(d['rank'] is not None and d['rank']<=n for d in regions)
        recall.append(dict(N=n,control_variant_hits=int(hits[:n].sum()),precision_proxy=float(hits[:n].mean()),
            regions_recalled=r,region_recall=r/8,unique_region_yield=r/n,threshold=float(table.iloc[n-1].total_score)))
    return dict(regions=regions,curve=recall,panel_control_prevalence=float(hits.mean()))


def drift_simulation(Y,meta,frozen,table,reps):
    """[inference] Separate Brownian/jump genotype experiment; no injected selection."""
    rng=np.random.default_rng(SEED); K=frozen['R'].shape[0]; mechanism_results={}
    # Brownian slope variance for each actual observed-date regression design.
    coefficients=[]; informations=[]
    for g in range(K):
        d=meta[meta.deme_id.eq(g)]; times=np.sort(d.t.unique()); x=d.t.to_numpy()-d.t.mean()
        if len(times)<2 or x@x<=0:
            continue
        weights=np.array([x[d.t.to_numpy()==tt].sum()/(x@x) for tt in times])
        covariance=np.minimum.outer(times-times[0],times-times[0])
        coefficients.append(float(weights@covariance@weights)); informations.append(float(x@x))
    weights=np.asarray(informations)/sum(informations)
    common_coefficient=float(np.sum(np.asarray(coefficients)*weights**2))
    diffusion=frozen['tau2']/common_coefficient
    chosen=rng.integers(0,Y.shape[1],size=reps)
    for mechanism in ['random_walk','random_walk_jump','pure_replacement']:
        S=np.empty((len(meta),reps)); baseline=[]
        for g in range(K):
            part=Y[meta.deme_id.eq(g)][:,chosen]
            baseline.append(np.nanmean(part,axis=0) if len(part) else np.full(reps,.5))
        baseline=np.nan_to_num(np.array(baseline),nan=.5)
        present=np.array([g for g in range(K) if meta.deme_id.eq(g).any()])
        donor=np.arange(K); donor[present]=rng.permutation(present)
        jumps=baseline[donor]-baseline
        sign=rng.choice([-1,1],size=reps)
        for g in range(K):
            ix=np.flatnonzero(meta.deme_id.to_numpy()==g); t=meta.t.to_numpy()[ix]
            if len(ix)==0: continue
            times=np.sort(np.unique(t)); path=np.repeat(baseline[g][None,:],len(times),axis=0)
            af=baseline[g]
            if mechanism!='pure_replacement':
                increments=rng.normal(size=(len(times)-1,reps))*np.sqrt(np.diff(times)[:,None]*diffusion*af*(1-af))
                path[1:]+=np.cumsum(increments,axis=0)
            if mechanism!='random_walk':
                when=rng.integers(1,len(times),size=reps)
                path+=(np.arange(len(times))[:,None]>=when)*abs(jumps[g])[None,:]*sign
            probabilities=np.clip(path[np.searchsorted(times,t)],.001,.999)
            S[ix]=rng.binomial(2,probabilities)/2
        result=score_matrix(S,meta,frozen)
        result['mechanism']=mechanism
        result.to_csv(OUT/f'simulation_{mechanism}.tsv',sep='\t',index=False)
        valid=result[result.valid]
        mechanism_results[mechanism]=dict(replicates=reps,valid=len(valid),diffusion_per_kyr=diffusion,
            score_quantiles={str(q):float(valid.total_score.quantile(q)) for q in [.5,.9,.95,.99,1]},
            shape_quantiles={str(q):float(valid['shape'].quantile(q)) for q in [.1,.5,.9]},
            rates=[dict(N=n,rate=float((valid.total_score>=table.iloc[n-1].total_score).mean()),
                exceedances=int((valid.total_score>=table.iloc[n-1].total_score).sum())) for n in ns_for(len(table))])
    return mechanism_results


def bootstrap(Y,meta,frozen,table,reps):
    rng=np.random.default_rng(SEED+1); original=table.index.to_numpy()
    ns=ns_for(len(table)); memberships={n:np.zeros(len(table),int) for n in ns}; metrics=[]
    strata=[d.index.to_numpy() for _,d in meta.groupby(['deme_id','node_id'])]
    base_rank=table.set_index(table.index)['rank']
    for b in range(reps):
        ix=np.concatenate([rng.choice(s,size=len(s),replace=True) for s in strata])
        result=score_matrix(Y[ix],meta.iloc[ix].reset_index(drop=True),frozen)
        ordered=result[result.valid].sort_values('total_score',ascending=False).index.to_numpy()
        row=dict(replicate=b,valid=len(ordered)); ranks=pd.Series(np.arange(1,len(ordered)+1),index=ordered)
        shared=base_rank.index.intersection(ranks.index)
        row['rank_spearman']=float(spearmanr(base_rank.loc[shared],ranks.loc[shared]).statistic)
        for n in ns:
            kept=np.isin(original,ordered[:n]); memberships[n]+=kept
            intersection=int(kept[:n].sum()); row[f'intersection_{n}']=intersection; row[f'jaccard_{n}']=intersection/(2*n-intersection)
        metrics.append(row)
        print(f'[observation] bootstrap {b+1}/{reps}',flush=True)
    stability=table[['locus','rank']].copy()
    for n in ns: stability[f'membership_top_{n}']=memberships[n]/reps
    stability.to_csv(OUT/'bootstrap_membership.tsv',sep='\t',index=False)
    pd.DataFrame(metrics).to_csv(OUT/'bootstrap_metrics.tsv',sep='\t',index=False)
    return dict(replicates=reps,mean_rank_spearman=float(pd.DataFrame(metrics).rank_spearman.mean()),
        metrics=[dict(N=n,mean_jaccard=float(np.mean([r[f'jaccard_{n}'] for r in metrics])),
            mean_intersection=float(np.mean([r[f'intersection_{n}'] for r in metrics])),
            mean_original_top_membership=float((memberships[n][:n]/reps).mean())) for n in ns])


def main():
    global OUT
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('--reuse',action='store_true')
    ap.add_argument('--bootstrap',type=int,default=30); ap.add_argument('--simulations',type=int,default=300)
    ap.add_argument('--panel-file',type=str,default=None,help='external variant panel TSV (locus,chrom,pos,akbari_X,pgen_index)')
    ap.add_argument('--out-dir',type=str,default=None,help='alternate output directory (for sharded runs)')
    args=ap.parse_args()
    if args.out_dir: OUT=Path(args.out_dir).resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    frozen=frozen_inputs()
    if args.reuse:
        panel=pd.read_csv(OUT/'panel.tsv',sep='\t'); meta=pd.read_csv(OUT/'sample_metadata.tsv',sep='\t'); Y=np.load(OUT/'genotypes.npy')
        inventory=json.loads((OUT/'inventory.json').read_text())
    else:
        panel,meta,Y,inventory=panel_and_data(args.panel_file); (OUT/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    inventory['verified_pvar_indices']=verify_variant_map(panel)
    (OUT/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print('[observation] scoring '+str(Y.shape),flush=True)
    scores=score_matrix(Y,meta,frozen); table=rank_table(panel,scores)
    table.to_csv(OUT/'priority_loci.tsv',sep='\t',index=False)
    control=validate_controls(table)
    result=dict(inventory=inventory,scored_variants=len(table),null='empirical common-slope variance law',
        frozen={k:v for k,v in frozen.items() if k!='R'},controls=control)
    result_path=OUT/'metrics.json'
    result_path.write_text(json.dumps(result,indent=2)+'\n')
    decile=int(np.ceil(len(table)/10)); top=table.iloc[:decile]; sig=table.akbari_abs_X.gt(5.45)
    novel=top[top.akbari_abs_X.le(5.45)]
    novel.to_csv(OUT/'top_decile_not_akbari_significant.tsv',sep='\t',index=False)
    result['akbari']=dict(covered_significant=int(sig.sum()),top_decile_significant=int(top.akbari_abs_X.gt(5.45).sum()),
        significant_in_top_decile_fraction=float(top.akbari_abs_X.gt(5.45).sum()/max(sig.sum(),1)),
        top_decile_not_significant=len(novel),spearman=float(spearmanr(table.total_score,table.akbari_abs_X).statistic),
        full_cache_significant=inventory['cache_significant'])
    result['carriage_sensitivity']={col:dict(spearman=float(spearmanr(table['rank'],table['rank_'+col]).statistic),
        top_decile_overlap=float((table.iloc[:decile]['rank_'+col]<=decile).mean())) for col in ['score_carriage_half','score_no_carriage']}
    result_path.write_text(json.dumps(result,indent=2)+'\n')
    print('[observation] controls '+json.dumps(control),flush=True)
    result['drift_simulation']=drift_simulation(Y,meta,frozen,table,args.simulations)
    result_path.write_text(json.dumps(result,indent=2)+'\n')
    if args.bootstrap:
        result['bootstrap']=bootstrap(Y,meta,frozen,table,args.bootstrap)
        result_path.write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4)); curve=control['curve']; ns=[r['N'] for r in curve]
    axes[0].plot(ns,[r['precision_proxy'] for r in curve],'o-',label='Control precision proxy')
    axes[0].plot(ns,[r['region_recall'] for r in curve],'o-',label='Region recall')
    axes[0].axhline(control['panel_control_prevalence'],ls=':',color='gray'); axes[0].legend()
    for mechanism,d in result['drift_simulation'].items(): axes[1].plot(ns,[r['rate'] for r in d['rates']],'o-',label=mechanism)
    for ax in axes: ax.set_xscale('log'); ax.set_xlabel('Top N variants')
    axes[1].set_ylabel('Simulated null exceeds Nth score'); axes[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(ROOT/'fig/prioritization_validation.png',dpi=160); plt.close(fig)
    sys.path.insert(0,str(ROOT)); import scripts.provenance as provenance
    provenance.get_git_info=lambda: {'disabled':'operator prohibits all git commands'}
    rec=provenance.RunRecorder('prioritization',command=['.venv/bin/python','analysis/prioritization_score.py',*sys.argv[1:]])
    rec.set_seed(SEED)
    for p in ['docs/prioritization-prereg.md','analysis/prioritization_score.py','analysis/mixture_null/genomewide_corrected_results.json']:
        rec.record_input(p)
    for p in ['data/derived/tier2_pool_cache.pkl','data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.pgen']:
        rec.manifest['inputs'][p]={'size_bytes':(ROOT/p).stat().st_size,'hash':'omitted: bounded immutable-input I/O'}
    for p in ['analysis/deme_rebuild/nodes.tsv','analysis/deme_rebuild/individuals.tsv','analysis/deme_rebuild/panels.npz',
              'analysis/mixture_null/genomewide_demes_akbari_pgen.tsv','data/derived/prioritization/panel.tsv']:
        rec.record_input(p)
    rec.record_parameter('frozen',result['frozen']); rec.record_parameter('bootstrap',args.bootstrap); rec.record_parameter('simulations',args.simulations)
    rec.record_parameter('verified_pvar_indices',inventory['verified_pvar_indices'])
    rec.record_output('data/derived/priority_loci.tsv'); rec.record_output('data/derived/prioritization/metrics.json')
    for p in OUT.glob('*.tsv'): rec.record_output(str(p.relative_to(ROOT)))
    rec.record_output('fig/prioritization_validation.png'); rec.save()
    print('[observation] complete '+str(len(table)),flush=True)

if __name__=='__main__': main()
