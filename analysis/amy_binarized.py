#!/usr/bin/env python3
"""Prespecified binary structural alleles; docs/amy-binarized-prereg.md."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import time
import zipfile
import numpy as np
import pandas as pd
from scipy import stats
from cn_as_snp import fit_row, empirical_p, ordering_gate, match_samples, AkbariGLMM
from load_akbari_geno import AkbariGeno, iter_pvar, paths, parse_psam
from load_akbari import iter_chunks
from mixture_null import fit, design

OUT = Path('data/derived/amy_binarized')
PRE = Path('docs/amy-binarized-prereg.md')
OLD = Path('data/derived/cn_as_snp')
SEED = 20261009
CUTS = {'primary_le1': 1, 'sensitivity_le2': 2, 'sensitivity_le3': 3}


def hap_copies(label):
    m = re.fullmatch(r'H(\d+)(?:\^[ar]|A\d+(?:B\d+)?)?', str(label))
    if not m:
        raise ValueError(f'unknown structural label: {label}')
    return int(m.group(1))


def dosage(frame, cut):
    return (frame.H1.map(hap_copies).gt(cut).astype(int) +
            frame.H2.map(hap_copies).gt(cut).astype(int)).to_numpy()


def write_json(name, obj):
    (OUT/name).write_text(json.dumps(obj, indent=2, default=str)+'\n')


def channels(age, y, groups):
    t = -np.asarray(age)/1000
    y = np.asarray(y)/2
    f = fit(t, y, groups)
    g, n, x = design(t, groups)
    means = np.bincount(g, weights=y)/n
    xp = t-t.mean()
    between = float(xp@means[g]/(xp@xp))
    pooled = float(xp@y/(xp@xp))
    contribution = float(f['beta']*(x@x)/(xp@xp))
    share = abs(between)/(abs(between)+abs(contribution))
    rng = np.random.default_rng(SEED+1)
    strata = [np.flatnonzero(g==j) for j in range(len(n))]
    bs = []
    for _ in range(1000):
        ix = np.concatenate([rng.choice(s,len(s),replace=True) for s in strata])
        tt, yy, gg = t[ix], y[ix], g[ix]
        mu = np.bincount(gg,weights=yy)/n
        xx = tt-tt.mean()
        bs.append(float(xx@mu[gg]/(xx@xx)))
    factor = stats.t.ppf(.975,f['df'])+stats.norm.ppf(.8)
    ideal_se = np.sqrt(np.sum(x*x*means[g]*(1-means[g])/2))/(x@x)
    return dict(beta_within=f['beta'], beta_between=between,
                beta_total=f['beta']+between, beta_pooled=pooled,
                within_pooled_contribution=contribution,
                within_HC3_se=f['se'], within_p=f['p'],
                between_bootstrap_se=float(np.std(bs,ddof=1)),
                carriage_absolute_share=share,
                channel_class='mostly carriage/composition' if share>=.75 else 'mostly within-history' if share<=.25 else 'both',
                opposite_signs=bool(between*contribution<0),
                n=len(y), B_frequency=float(y.mean()), df=f['df'],
                within_time_information=float(x@x),
                MDE80_HC3_per_kyr=float(factor*f['se']),
                MDE80_ideal_binomial_per_kyr=float(factor*ideal_se),
                MDE80_HC3_2kyr=float(2*factor*f['se']),
                MDE80_HC3_full_span=float(np.ptp(t)*factor*f['se']),
                date_span_kyr=float(np.ptp(t)),
                observed_within_2kyr=float(2*f['beta']))


def freq_rows(frame, scope, cut_name, cut, dimension='component'):
    y = dosage(frame,cut)
    labels = frame.p2 if dimension=='component' else pd.cut(frame.age,[0,2000,4000,6000,8000,10000,np.inf],right=False)
    result=[]
    for label in labels.dropna().unique():
        ix=(labels==label).to_numpy(); yy=y[ix]
        result.append(dict(scope=scope,cut=cut_name,dimension=dimension,group=str(label),n=int(ix.sum()),
                           mean_BP=float(frame.loc[ix,'age'].mean()),B_copies=int(yy.sum()),alleles=int(2*len(yy)),
                           B_frequency=float(yy.mean()/2),n0=int((yy==0).sum()),n1=int((yy==1).sum()),n2=int((yy==2).sum())))
    return sorted(result,key=lambda r:-r['mean_BP'])


def extract():
    scratch=OUT/'scratch'; scratch.mkdir(parents=True,exist_ok=True)
    tables={}
    with zipfile.ZipFile('data/raw/amy/amylase_diversity_project-v1.1.zip') as z:
        for suffix in ['haplotype_structure_for_supplement.tsv','ancient_metadata_for_supplement.tsv','copy_number_table_for_supplement.tsv']:
            names=[n for n in z.namelist() if n.endswith('/'+suffix)]
            if len(names)!=1: raise ValueError(suffix)
            p=scratch/suffix; p.write_bytes(z.read(names[0])); tables[suffix]=pd.read_csv(p,sep='\t')
    return tables


def copy_check(h, cn):
    merged=h.merge(cn[['sample','AMY1']],on='sample',validate='one_to_one')
    merged=merged.dropna(subset=['H1','H2','AMY1'])
    labels=sorted(set(merged.H1)|set(merged.H2))
    X=np.column_stack([(merged.H1.eq(l).astype(int)+merged.H2.eq(l).astype(int)).to_numpy() for l in labels])
    estimates,_,rank,_=np.linalg.lstsq(X,merged.AMY1.to_numpy(float),rcond=None)
    rows=[]
    for label,est in zip(labels,estimates):
        homo=merged[merged.H1.eq(label)&merged.H2.eq(label)]
        rows.append(dict(haplotype=label,name_AMY1_copies=hap_copies(label),homozygote_n=len(homo),
                         homozygote_mean_AMY1_per_hap=float(homo.AMY1.mean()/2) if len(homo) else None,
                         additive_estimated_AMY1_per_hap=float(est),diplotype_slots=int(X[:,labels.index(label)].sum())))
    pd.DataFrame(rows).to_csv(OUT/'haplotype_copy_check.tsv',sep='\t',index=False)
    return dict(n=len(merged),rank=int(rank),labels=len(labels),RMSE=float(np.sqrt(np.mean((merged.AMY1-X@estimates)**2))))


def provenance(mode, start):
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    from scripts import provenance as p
    p.get_git_info=lambda: {'commit':'not inspected: no-git rule','branch':'not inspected','dirty':'unknown'}
    r=p.RunRecorder('amy-binarized-'+mode,command=[sys.executable]+sys.argv)
    r.set_seed(SEED); r.record_parameter('start_utc',start)
    for path in [PRE,Path(__file__),OLD/'primary_cohort.tsv',OLD/'summary.json',OLD/'GRM_LOCO_chr1.npy',Path('sim/akbari_glmm.py'),Path('data/raw/amy/amylase_diversity_project-v1.1.zip')]: r.record_input(path)
    for path in paths()[:2]+(Path('data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz'),):
        st=path.stat(); r.record_parameter(str(path),{'size':st.st_size,'mtime_ns':st.st_mtime_ns,'sha256':'not computed: bounded run'})
    for path in OUT.glob('*'):
        if path.is_file(): r.record_output(path)
    r.save()


def biological(seconds=1000):
    start=datetime.now(timezone.utc).isoformat(); deadline=time.monotonic()+seconds
    state=dict(started_utc=start,prereg_mtime_ns=PRE.stat().st_mtime_ns,prereg_sha256=hashlib.sha256(PRE.read_bytes()).hexdigest(),status='started',results=[],nulls=[])
    write_json('summary.json',state)
    tables=extract(); h=tables['haplotype_structure_for_supplement.tsv']; cn=tables['copy_number_table_for_supplement.tsv']; meta=tables['ancient_metadata_for_supplement.tsv']
    state['copy_check']=copy_check(h,cn)
    ancient=h[h.type.eq('Ancient')].merge(meta[['sample','age']],on='sample',validate='one_to_one')
    ancient['age']=pd.to_numeric(ancient.age,errors='coerce'); ancient=ancient.dropna(subset=['age','H1','H2','p2'])
    cohort=pd.read_csv(OLD/'primary_cohort.tsv',sep='\t')
    mapped=match_samples(ancient,parse_psam(paths()[2])); matched=mapped[mapped.pgen_full_index.notna()].sort_values('pgen_full_index').reset_index(drop=True)
    assert list(matched['sample'])==list(cohort['sample'])
    assert np.array_equal(matched.age,cohort.age)
    for col in ['H1','H2','p2']: assert list(matched[col])==list(cohort[col])
    state['cohort']=dict(labellable=len(ancient),matched=len(cohort),absent=int(mapped.pgen_full_index.isna().sum()))
    cohort.to_csv(OUT/'primary_cohort.tsv',sep='\t',index=False)
    descriptions=[]; phen={}; state['channels']={}; counts={}
    for name,cut in CUTS.items():
        y=dosage(cohort,cut); phen[name]=y; counts[name]=dict(n=len(y),n0=int((y==0).sum()),n1=int((y==1).sum()),n2=int((y==2).sum()),B_copies=int(y.sum()),alleles=2*len(y),B_frequency=float(y.mean()/2))
        state['channels'][name]=channels(cohort.age,y,cohort.p2)
        for scope,df in [('matched_225',cohort),('labelled_288',ancient)]:
            for dim in ['component','date_bin']: descriptions+=freq_rows(df,scope,name,cut,dim)
    state['dosage_counts']=counts
    checks=[]
    for label,part in ancient.groupby('p2'):
        checks.append(dict(component=label,n=len(part),mean_BP=float(part.age.mean()),
                           H1a_frequency=float((part.H1.eq('H1^a').sum()+part.H2.eq('H1^a').sum())/(2*len(part))),
                           H5_frequency=float((part.H1.eq('H5').sum()+part.H2.eq('H5').sum())/(2*len(part)))))
    pd.DataFrame(checks).sort_values('mean_BP',ascending=False).to_csv(OUT/'verified_ancient_haplotypes.tsv',sep='\t',index=False)
    pd.DataFrame(descriptions).to_csv(OUT/'frequencies.tsv',sep='\t',index=False)
    pd.DataFrame({'sample':cohort['sample'],**phen}).to_csv(OUT/'dosages.tsv',sep='\t',index=False)
    modern=h[h.type.eq('Modern')].dropna(subset=['H1','H2','p1']); modernrows=[]
    for name,cut in CUTS.items():
        for label,part in modern.groupby('p1'):
            y=dosage(part,cut)
            modernrows.append(dict(cut=name,population=label,n=len(y),B_copies=int(y.sum()),alleles=2*len(y),B_frequency=float(y.mean()/2),n0=int((y==0).sum()),n1=int((y==1).sum()),n2=int((y==2).sum()),
                                   H5_frequency=float((part.H1.eq('H5').sum()+part.H2.eq('H5').sum())/(2*len(part))),
                                   H9_frequency=float((part.H1.eq('H9').sum()+part.H2.eq('H9').sum())/(2*len(part)))))
    pd.DataFrame(modernrows).to_csv(OUT/'modern_frequencies.tsv',sep='\t',index=False)
    # Recheck carried-forward CN means and within slopes on full measured ancients.
    measured=cn[cn.type.eq('Ancient')].merge(meta[['sample','age']],on='sample',validate='one_to_one')
    measured['age']=pd.to_numeric(measured.age,errors='coerce'); measured=measured.dropna(subset=['age','AMY1','p2'])
    verification=[]
    for label,part in measured.groupby('p2'):
        verification.append(dict(component=label,mean_AMY1=float(part.AMY1.mean()),**fit(-part.age/1000,part.AMY1,part.p2)))
    pd.DataFrame(verification).to_csv(OUT/'verified_CN.tsv',sep='\t',index=False)
    (Path('docs/amy-binarized-results.md')).open('a').write('\n## NEW INCREMENT 1 — '+start+'\n\n[observation] Preregistration persisted before this run. Reconstructed '+str(len(ancient))+' labelled ancients; exact prior cohort '+str(len(cohort))+' verified. Frequency, copy-check, modern, and power/channel artifacts now saved; fitted results pending.\n')
    write_json('summary.json',state)
    previous=json.loads((OLD/'summary.json').read_text()); refs=previous['references']; positions={}; needed={r['ID'] for r in refs.values()}
    for i,(_,_,vid,_,_) in enumerate(iter_pvar(paths()[1])):
        if vid in needed: positions[vid]=i
        if len(positions)==len(needed): break
    Ks={c:np.load(OLD/f'GRM_LOCO_chr{c}.npy') for c in [1,2,5]}; zero=np.zeros_like(Ks[1]); dates=cohort.age.to_numpy(float)
    controls=[]
    with AkbariGeno(sample_subset=cohort.pgen_full_index.to_numpy(dtype=np.uint32)) as geno:
        bg=pd.read_csv(OLD/'background_markers.tsv',sep='\t')
        genotypes=geno.read_variants(bg['index'].to_numpy(dtype=np.uint32)).T.astype(float)
        genotypes[(genotypes<0)|(genotypes>2)]=np.nan
        for chrom in Ks:
            rebuilt=AkbariGLMM.compute_grm(genotypes,marker_chromosomes=bg.chrom.to_numpy(),target_chromosome=chrom)
            np.testing.assert_allclose(rebuilt,Ks[chrom],rtol=1e-12,atol=1e-12)
        state['GRM_verified']=dict(background_markers=len(bg),AMY_LOCO_markers=int(bg.chrom.ne(1).sum()),chromosomes=int(bg.chrom.nunique()),rebuilt_chromosomes=list(Ks))
        for label,ref in refs.items():
            yy=geno.read_hardcalls(np.array([positions[ref['ID']]],dtype=np.uint32))[0].astype(float); yy[yy<0]=np.nan
            for model,k in [('GRM_on',Ks[int(ref['CHROM'])]),('GRM_off',zero)]:
                r=fit_row(label,model,yy,dates,k); r['released_X']=ref['X']; controls.append(r)
    state['controls']=controls; state['ordering_gate_passed']=ordering_gate(controls)
    write_json('summary.json',state)
    if not state['ordering_gate_passed']: raise RuntimeError('previous ordering gate not reproduced')
    for name,y in phen.items():
        for model,k in [('GRM_on',Ks[1]),('GRM_off',zero)]:
            r=fit_row(name,model,y,dates,k)
            if r.get('converged'):
                for anchor in ['LCT','SLC45A2']:
                    a=next(q for q in controls if q['locus']==anchor and q['model']==model)
                    r[anchor+'_extrapolated_abs_X']=abs(a['released_X'])*abs(r['z_score'])/abs(a['z_score'])
            state['results'].append(r)
    write_json('summary.json',state); pd.DataFrame(state['results']).to_csv(OUT/'statistics.tsv',sep='\t',index=False)
    print(json.dumps(state['results'],indent=2),flush=True)
    rng=np.random.default_rng(SEED); strata=[np.flatnonzero(cohort.p2.eq(g)) for g in cohort.p2.unique()]
    observed={(r['locus'],r['model']):r for r in state['results'] if r.get('converged')}
    null={}; failures={}
    with (OUT/'permutations.tsv').open('w') as f:
        f.write('null\tcut\tmodel\tpermutation\tz_score\tconverged\terror\n')
        for nulltype in ['ancestry','unrestricted']:
            for b in range(199):
                if time.monotonic()>deadline: break
                for name,y in phen.items():
                    perm=y.copy()
                    for s in strata if nulltype=='ancestry' else [np.arange(len(y))]: perm[s]=y[rng.permutation(s)]
                    for model,k in [('GRM_on',Ks[1]),('GRM_off',zero)]:
                        if time.monotonic()>deadline or (name,model) not in observed: continue
                        key=(nulltype,name,model); values=null.setdefault(key,[]); failures.setdefault(key,0)
                        r=fit_row(name,model,perm,dates,k)
                        if r.get('converged'): values.append(r['z_score'])
                        else: failures[key]+=1
                        f.write(f"{nulltype}\t{name}\t{model}\t{b}\t{r.get('z_score','')}\t{r.get('converged',False)}\t{r.get('error','')}\n"); f.flush()
                state['nulls']=[dict(null=key[0],cut=key[1],model=key[2],requested=199,successful=len(vals),failed=failures[key],empirical_P=empirical_p(observed[(key[1],key[2])]['z_score'],vals),null_abs_Z_quantiles=np.quantile(abs(np.asarray(vals)),[.5,.9,.95,.99]).tolist() if vals else []) for key,vals in null.items()]
                write_json('summary.json',state)
                if (b+1)%20==0: print(nulltype,b+1,flush=True)
    state['status']='complete' if len(state['nulls'])==12 and all(r['successful']+r['failed']==199 for r in state['nulls']) else 'partial'
    state['finished_utc']=datetime.now(timezone.utc).isoformat(); write_json('summary.json',state)
    provenance('biology',start)


def distribution():
    start=datetime.now(timezone.utc).isoformat()
    exact=json.loads(Path('data/derived/two_channel/akbari_audit.json').read_text())
    thresholds={'LCT':16.523199593858344,'FADS':16.153410398,'SLC45A2':13.149626229714274,'SLC24A5':5.990305093,'AMY':3.794956889120}
    # Exact maxima from independently streamed previous window audit, not rounded substitutes.
    for name in thresholds:
        key='AMY1' if name=='AMY' else name
        thresholds[name]=exact[key]['max_abs_X']
    thresholds.update({f'rounded_{k}':v for k,v in zip(['LCT','FADS','SLC45A2','SLC24A5','AMY'],[16.52,16.15,13.15,5.99,3.795])})
    windows={'AMY':(1,103500000,104500000),'LCT':(2,136300000,136900000),'FADS':(11,61300000,61700000),'SLC45A2':(5,33800000,34200000),'SLC24A5':(15,48200000,48600000)}
    audit={k:dict(n=0,max_abs_X=0.,strong_n=0) for k in windows}
    counts={k:0 for k in thresholds}; parts=[]; bandblocks=set(); strongblocks=set(); maxima={}; missing=0; band=0; strong=0; total=0; top=[]
    for d in iter_chunks(usecols=['CHROM','POS','ID','X']):
        ok=np.isfinite(d.X); missing+=int((~ok).sum()); d=d[ok]; a=d.X.abs().to_numpy(); parts.append(a); total+=len(a)
        for name,(chrom,lo,hi) in windows.items():
            region=d[d.CHROM.eq(chrom)&d.POS.between(lo,hi)]
            if len(region):
                audit[name]['n']+=len(region)
                audit[name]['strong_n']+=int(region.X.abs().gt(5.45).sum())
                audit[name]['max_abs_X']=max(audit[name]['max_abs_X'],float(region.X.abs().max()))
        top.extend(d.loc[d.X.abs().nlargest(20).index,['ID','X']].to_dict('records'))
        top=sorted(top,key=lambda r:-abs(r['X']))[:20]
        bm=(a>=3.16)&(a<=5.45); sm=a>5.45; band+=int(bm.sum()); strong+=int(sm.sum())
        for k,v in thresholds.items(): counts[k]+=int((a>v).sum())
        for target,mask in [(bandblocks,bm),(strongblocks,sm)]:
            target.update(zip(d.loc[mask,'CHROM'].astype(int), (d.loc[mask,'POS']//250000).astype(int)))
        blockmax=d.assign(block=d.POS//250000,abs_X=a).groupby(['CHROM','block']).abs_X.max()
        for key,v in blockmax.items(): maxima[key]=max(v,maxima.get(key,0))
        write_json('distribution_partial.json',dict(n_finite=total,n_band=band,n_strong=strong))
    a=np.concatenate(parts)
    output=dict(started_utc=start,n_finite=total,n_missing=missing,quantile_probabilities=[0,.25,.5,.75,.9,.95,.99,.999,1],quantiles=np.quantile(a,[0,.25,.5,.75,.9,.95,.99,.999,1]).tolist(),thresholds={k:dict(abs_X=v,n_exceed=counts[k],fraction_exceed=counts[k]/total) for k,v in thresholds.items()},band_variants=band,strong_variants=strong,band_any_blocks=len(bandblocks),strong_any_blocks=len(strongblocks),overlap_blocks=len(bandblocks&strongblocks),band_only_any_blocks=len(bandblocks-strongblocks),all_blocks=len(maxima),band_max_blocks=sum(3.16<=v<=5.45 for v in maxima.values()),strong_max_blocks=sum(v>5.45 for v in maxima.values()))
    output['top_variants']=top
    output['window_audit']=audit
    for name in windows: assert audit[name]['max_abs_X']==thresholds[name]
    assert audit['AMY']['n']==2204 and audit['AMY']['strong_n']==0
    pd.DataFrame([dict(chrom=c,block=int(b),max_abs_X=v) for (c,b),v in maxima.items()]).to_csv(OUT/'genomewide_block_maxima.tsv',sep='\t',index=False)
    write_json('distribution.json',output); print(json.dumps(output,indent=2),flush=True); provenance('distribution',start)


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--mode',choices=['biology','distribution'],required=True); p.add_argument('--seconds',type=int,default=1000); args=p.parse_args()
    if not PRE.exists(): raise RuntimeError('write preregistration before results')
    OUT.mkdir(parents=True,exist_ok=True)
    biological(args.seconds) if args.mode=='biology' else distribution()
