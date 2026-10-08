#!/usr/bin/env python3
"""Bounded genome-wide-distributed measured-label residual control scan.

Protocol: docs/mixture-null-prereg.md. No biological interpretation before
held-out calibration. This is not fractional-ancestry correction.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from load_aadr import load_bundle, iter_snp, MAIN_PREFIX, HO_PREFIX, ANNO_NAME, HO_ANNO_NAME
from load_akbari_geno import AkbariGeno, iter_pvar, paths
from load_akbari import iter_chunks
from mixture_null import fit

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/mixture_null'


def metadata(prefix, anno_name):
    gm,ind,anno=load_bundle(prefix=prefix,anno_name=anno_name)
    seen=set(); records=[]
    for i,r in enumerate(ind):
        a=anno.get(r.id)
        if a is None or a.is_modern or a.date_bp is None or not 0<a.date_bp<=14000:
            continue
        if any(word in a.assessment.upper() for word in ['FAIL','CRITICAL']):
            continue
        if a.individual_id in seen:
            continue
        seen.add(a.individual_id)
        records.append(dict(row=i,id=r.id,individual=a.individual_id,group=r.group,t=-a.date_bp/1000))
    df=pd.DataFrame(records)
    groups=df.groupby('group').agg(n=('t','size'),distinct_dates=('t','nunique'),lo=('t','min'),hi=('t','max'))
    groups['span_kyr']=groups.hi-groups.lo
    groups['eligible']=(groups.n>=16)&(groups.distinct_dates>=3)&(groups.span_kyr>=2)
    return gm,df[df.group.isin(groups.index[groups.eligible])].copy(),groups


def eligible_snps(prefix):
    result={}
    for i,s in enumerate(iter_snp(str(ROOT/'data/raw/aadr-v66.p1'/f'{prefix}.snp'))):
        if i%200!=0 or s[1] in ['1','6'] or s[1] not in [str(j) for j in range(1,23)]:
            continue
        key=(s[1],int(s[3]),s[4].upper(),s[5].upper())
        result[key]=(i,s[0])
    return result


def matched_sites(eligible, limit):
    matches=[]
    _,pvar,_=paths()
    for pi,(chrom,pos,vid,ref,alt) in enumerate(iter_pvar(pvar)):
        key=(chrom,pos,ref.upper(),alt.upper())
        if key in eligible:
            ai,rsid=eligible[key]
            matches.append(dict(pgen_index=pi,aadr_index=ai,id=vid,rsid=rsid,chrom=chrom))
    if len(matches)>limit:
        matches=[matches[i] for i in np.linspace(0,len(matches)-1,limit,dtype=int)]
    return pd.DataFrame(matches)


def score_matrix(Y,meta,sites,label):
    results=[]
    t=meta.t.to_numpy(); groups=meta.group.to_numpy()
    for j,row in sites.reset_index(drop=True).iterrows():
        y=np.asarray(Y[:,j],float)
        good=np.isfinite(y)
        supported=[]
        for g in np.unique(groups[good]):
            w=good&(groups==g)
            if np.sum(w)>=8 and len(np.unique(t[w]))>=3 and np.ptp(t[w])>=2:
                supported.append(g)
        ok=good&np.isin(groups,supported)
        if np.sum(ok)<100 or len(supported)<2:
            continue
        af=float(y[ok].mean())
        if not .1<=af<=.9:
            continue
        f=fit(t[ok],y[ok],groups[ok])
        if not f['identified'] or not np.isfinite(f['p']):
            continue
        results.append(dict(bundle=label,id=row.id,chrom=row.chrom,af=af,**f))
    return pd.DataFrame(results)


def report(df,label):
    if len(df)<2:
        return dict(bundle=label,tested=len(df),calibration_gate=False,reason='insufficient controls')
    df=df.copy()
    df['partition']=np.where(np.arange(len(df))<len(df)//2,'train','heldout')
    train=df[df.partition.eq('train')]; held=df[df.partition.eq('heldout')]
    count=int((held.p<=.05).sum()); N=len(held)
    ci=stats.binomtest(count,N).proportion_ci()
    result=dict(bundle=label,tested=len(df),heldout_n=N,heldout_rejects=count,
                control_pass_fraction=count/N,ci_low=float(ci.low),ci_high=float(ci.high),
                calibration_gate=N>=100 and ci.low<=.05<=ci.high,
                floor_unit='allele-frequency/kyr',
                raw_mde80_quantiles={str(q):float(df.mde80.quantile(q)) for q in [.1,.5,.9]},
                chromosomes=sorted(df.chrom.astype(int).unique().tolist()),
                median_n=float(df.n.median()),median_k=float(df.k.median()))
    df.to_csv(OUT/f'genomewide_{label}.tsv',sep='\t',index=False)
    return result


def scan(limit=1200):
    OUT.mkdir(exist_ok=True)
    gm,meta,groups=metadata(MAIN_PREFIX,ANNO_NAME)
    groups.to_csv(OUT/'aadr1240k_groups.tsv',sep='\t')
    hgm,hmeta,hgroups=metadata(HO_PREFIX,HO_ANNO_NAME)
    hgroups.to_csv(OUT/'aadrHO_groups.tsv',sep='\t')
    inventory=dict(aadr1240k_groups=len(groups),aadr1240k_qualifying_groups=int(groups.eligible.sum()),
                   aadr1240k_qualifying_people=len(meta),aadrHO_groups=len(hgroups),
                   aadrHO_qualifying_groups=int(hgroups.eligible.sum()),aadrHO_qualifying_people=len(hmeta))
    (OUT/'genomewide_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print(json.dumps(inventory),flush=True)
    sites=matched_sites(eligible_snps(MAIN_PREFIX),limit)
    sites.to_csv(OUT/'genomewide_sites.tsv',sep='\t',index=False)
    print(f'exact matched sites {len(sites)}',flush=True)
    Y=np.array([gm.individual_genotypes(int(row),sites.aadr_index.to_numpy(int)) for row in meta.row],dtype=float)/2
    Y[Y>1]=np.nan
    results=[report(score_matrix(Y,meta,sites,'aadr1240k'),'aadr1240k')]
    print(json.dumps(results[-1]),flush=True)
    geno=AkbariGeno()
    sample_of={str(s):i for i,s in enumerate(geno.sample_ids)}
    pmeta=meta[meta.id.isin(sample_of)].copy()
    pmeta['pgen_row']=[sample_of[s] for s in pmeta.id]
    pmeta=pmeta.sort_values('pgen_row')
    geno.change_sample_subset(pmeta.pgen_row.to_numpy(np.uint32))
    P=geno.read_variants(sites.pgen_index.to_numpy(np.uint32)).T.astype(float)/2
    P[(P<0)|(P>1)]=np.nan
    geno.close()
    results.append(report(score_matrix(P,pmeta,sites,'akbari_pgen'),'akbari_pgen'))
    print(json.dumps(results[-1]),flush=True)
    hsites=matched_sites(eligible_snps(HO_PREFIX),limit)
    hsites.to_csv(OUT/'genomewide_HO_sites.tsv',sep='\t',index=False)
    H=np.array([hgm.individual_genotypes(int(row),hsites.aadr_index.to_numpy(int)) for row in hmeta.row],dtype=float)/2
    H[H>1]=np.nan
    results.append(report(score_matrix(H,hmeta,hsites,'aadrHO'),'aadrHO'))
    print(json.dumps(results[-1]),flush=True)
    summary={**inventory,'requested_controls':limit,'pgen_shared_people':len(pmeta),'scans':results}
    (OUT/'genomewide_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def deme_scan():
    """Secondary preregistered fallback; independently validate saved membership."""
    from load_aadr import parse_ind, parse_anno
    base=ROOT/'analysis/deme_rebuild'
    members=pd.read_csv(base/'individuals.tsv',sep='\t')
    paths_df=pd.read_csv(base/'paths.tsv',sep='\t')
    nodes=pd.read_csv(base/'nodes.tsv',sep='\t')
    gm,ind,anno=load_bundle()
    if members.biological_id.duplicated().any() or members.individual.duplicated().any():
        raise ValueError('deme identities duplicated')
    for r in members.itertuples():
        if ind[r.row].id != r.individual:
            raise ValueError('saved deme row ID mismatch')
        a=anno[r.individual]
        if a.individual_id!=r.biological_id or abs(a.date_bp-r.date_BP)>1e-6:
            raise ValueError('saved deme metadata mismatch')
    qualifying=paths_df[paths_df.qualifies].copy()
    if not (qualifying.n_nodes.ge(3)&qualifying.n.ge(30)&qualifying.span.ge(2000)).all():
        raise ValueError('path fails registered duration/count criteria')
    selected=members[members.deme_id.ge(0)].copy()
    for d,part in selected.groupby('deme_id'):
        node_part=nodes[nodes.deme_id.eq(d)]
        path=qualifying[qualifying.deme_id.eq(d)].iloc[0]
        bins=np.sort(node_part.time_bin.unique())
        if not np.all(np.diff(bins)==1) or len(bins)!=path.n_nodes:
            raise ValueError('nonadjacent or duplicated temporal nodes')
        if len(part)!=path.n or part.node_id.nunique()!=path.n_nodes or np.ptp(part.date_BP)<2000:
            raise ValueError('deme count/time mismatch')
        if part.lat_cell.nunique()!=1 or part.lon_cell.nunique()!=1:
            raise ValueError('deme crosses geographic cells')
        for node,p in part.groupby('node_id'):
            r=node_part[node_part.node_id.eq(node)].iloc[0]
            if len(p)!=r.n or abs(p.date_BP.mean()-r.date_BP)>1e-6:
                raise ValueError('node summary mismatch')
    selected['group']='deme_'+selected.deme_id.astype(int).astype(str)
    selected['id']=selected.individual
    selected['t']=-selected.date_BP/1000
    panel=np.load(base/'panels.npz')['embedding_sites']
    sites=pd.read_csv(OUT/'genomewide_sites.tsv',sep='\t')
    sites=sites[~sites.aadr_index.isin(panel)].copy()
    inventory=dict(qualifying_demes=len(qualifying),people=len(selected),nodes=len(nodes),
                   span_kyr_min=float(qualifying.span.min()/1000),span_kyr_max=float(qualifying.span.max()/1000),
                   n_min=int(qualifying.n.min()),n_max=int(qualifying.n.max()),
                   removed_embedding_controls=1200-len(sites),
                   validation='row IDs, dates, biological uniqueness, node counts, adjacent bins and geography verified; source path rules unit-tested')
    (OUT/'deme_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n')
    print(json.dumps(inventory),flush=True)
    Y=np.array([gm.individual_genotypes(int(row),sites.aadr_index.to_numpy(int)) for row in selected.row],float)/2
    Y[Y>1]=np.nan
    results=[report(score_matrix(Y,selected,sites,'demes_aadr1240k'),'demes_aadr1240k')]
    print(json.dumps(results[-1]),flush=True)
    geno=AkbariGeno(); sample_of={str(s):i for i,s in enumerate(geno.sample_ids)}
    pmeta=selected[selected.id.isin(sample_of)].copy()
    pmeta['pgen_row']=[sample_of[s] for s in pmeta.id]
    pmeta=pmeta.sort_values('pgen_row')
    geno.change_sample_subset(pmeta.pgen_row.to_numpy(np.uint32))
    P=geno.read_variants(sites.pgen_index.to_numpy(np.uint32)).T.astype(float)/2
    P[(P<0)|(P>1)]=np.nan; geno.close()
    results.append(report(score_matrix(P,pmeta,sites,'demes_akbari_pgen'),'demes_akbari_pgen'))
    print(json.dumps(results[-1]),flush=True)
    summary={**inventory,'pgen_shared_people':len(pmeta),'scans':results}
    (OUT/'genomewide_deme_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def calibrated_deme_controls():
    """One global temporal variance per bundle, fitted only on TRAIN loci."""
    result=[]
    rng=np.random.default_rng(20261008)
    for label in ['demes_aadr1240k','demes_akbari_pgen']:
        df=pd.read_csv(OUT/f'genomewide_{label}.tsv',sep='\t')
        train=df[df.partition.eq('train')]; held=df[df.partition.eq('heldout')].copy()
        moments=((train.beta.to_numpy()**2-train.se.to_numpy()**2)/
                 (train.af.to_numpy()*(1-train.af.to_numpy())))
        tau2=max(0.,float(moments.mean()))
        boots=np.maximum(0,rng.choice(moments,size=(500,len(moments)),replace=True).mean(axis=1))
        df['corrected_se']=np.sqrt(df.se**2+tau2*df.af*(1-df.af))
        df['corrected_p']=2*stats.t.sf(abs(df.beta/df.corrected_se),len(train)-1)
        cutoff=float(stats.t.ppf(.975,len(train)-1))
        df['corrected_mde80']=(cutoff+stats.norm.ppf(.8))*df.corrected_se
        held=df[df.partition.eq('heldout')]
        count=int((held.corrected_p<=.05).sum()); N=len(held)
        ci=stats.binomtest(count,N).proportion_ci()
        noise=rng.normal(size=(10000,N))
        powers=[]
        for effect in [0,.01,.02,.03,.04,.05,.075,.1]:
            powers.append(dict(effect=effect,unit='allele-frequency/kyr',
                               mean_conditional_power=float(np.mean(np.abs(noise+effect/held.corrected_se.to_numpy()[None,:])>cutoff))))
        r=dict(bundle=label,shared_demographic_variance_parameters=1,
               training_loci=len(train),tau2=tau2,
               tau2_bootstrap_quantiles={str(q):float(np.quantile(boots,q)) for q in [.025,.5,.975]},
               heldout_n=N,heldout_rejects=count,control_pass_fraction=count/N,
               ci_low=float(ci.low),ci_high=float(ci.high),calibration_gate=N>=100 and ci.low<=.05<=ci.high,
               floor_unit='allele-frequency/kyr',
               corrected_mde80_quantiles={str(q):float(held.corrected_mde80.quantile(q)) for q in [.1,.5,.9]},
               power=powers,seed=20261008,score_simulation_replicates=10000)
        df.to_csv(OUT/f'genomewide_corrected_{label}.tsv',sep='\t',index=False)
        result.append(r)
    (OUT/'genomewide_corrected_results.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def join_stats():
    # Joins raw summary release, not cached author-selected positive controls.
    ids=set()
    for p in OUT.glob('genomewide_*.tsv'):
        if p.stem in ['genomewide_aadr1240k','genomewide_akbari_pgen','genomewide_aadrHO','genomewide_demes_aadr1240k','genomewide_demes_akbari_pgen']:
            ids.update(pd.read_csv(p,sep='\t',usecols=['id']).id)
    parts=[]
    for chunk in iter_chunks(usecols=['ID','X','POSTERIOR','FDR','FILTER'],chunksize=500000):
        d=chunk[chunk.ID.isin(ids)]
        if len(d): parts.append(d)
    matches=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=['ID','X','POSTERIOR','FDR','FILTER'])
    matches.to_csv(OUT/'genomewide_akbari_stats.tsv',sep='\t',index=False)
    result=dict(unique_control_ids=len(ids),summary_matches=len(matches),
                median_X=float(matches.X.median()),abs_X_gt_5_45=int((matches.X.abs()>5.45).sum()),
                abs_X_ge_3_16=int((matches.X.abs()>=3.16).sum()),
                X_centered=False)
    (OUT/'genomewide_stats_join.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--controls',type=int,default=1200)
    ap.add_argument('--join-stats',action='store_true')
    ap.add_argument('--demes',action='store_true',help='registered validated longitudinal fallback')
    ap.add_argument('--calibrate-demes',action='store_true',help='TRAIN-only temporal variance repair')
    a=ap.parse_args()
    if a.join_stats: join_stats()
    elif a.calibrate_demes: print(json.dumps(calibrated_deme_controls(),indent=2))
    elif a.demes: print(json.dumps(deme_scan(),indent=2))
    else: print(json.dumps(scan(a.controls),indent=2))

if __name__=='__main__': main()
