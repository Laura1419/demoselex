#!/usr/bin/env python3
"""Measured-label (one-hot) mixture residual; protocol docs/mixture-null-prereg.md.

This is not fractional-ancestry correction or a CN evolutionary-process model.
All effects are forward copies/kyr; no unseen donor is fitted.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis' / 'mixture_null'
SEED = 20261008
EFFECTS = [0, .05, .1, .2, .3, .5, .75, 1, 1.5, 2]


def design(t, groups):
    t = np.asarray(t, float)
    _, g = np.unique(np.asarray(groups, str), return_inverse=True)
    n = np.bincount(g)
    x = t - (np.bincount(g, weights=t) / n)[g]
    return g, n, x


def fit(t, y, groups):
    """Measured means cancel exactly in slope; HC3 charges all group baselines."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    groups = np.asarray(groups, str)
    ok = np.isfinite(t) & np.isfinite(y)
    t, y, groups = t[ok], y[ok], groups[ok]
    g, n, x = design(t, groups)
    k, N, sxx = len(n), len(y), float(x @ x)
    result = dict(n=N, k=k, df=N-k-1, time_information=sxx)
    if sxx < 1e-12 or N-k-1 <= 0:
        return {**result, 'beta': None, 'se': None, 'p': None, 'mde80': None,
                'identified': False}
    mu = np.bincount(g, weights=y) / n
    beta = float(x @ (y-mu[g]) / sxx)
    e = y-mu[g]-beta*x
    h = 1/n[g] + x*x/sxx
    if np.any(h >= 1-1e-12):
        # Singleton observations contain no within-group information; their
        # zero residual cannot be divided by zero in HC3.
        useful = h < 1-1e-12
    else:
        useful = np.ones(N, bool)
    se = float(np.sqrt(np.sum((x[useful]*e[useful]/(1-h[useful]))**2))/sxx)
    se_iid = float(np.sqrt(e @ e / (N-k-1) / sxx))
    z = beta/se if se > 0 else np.nan
    p = float(2*stats.t.sf(abs(z), N-k-1))
    return {**result, 'beta': beta, 'se': se, 'se_iid': se_iid,
            'z': float(z), 'p': p, 'mde80': float((1.96+stats.norm.ppf(.8))*se),
            'identified': True}


def load_amy():
    raw = ROOT / 'data/raw/amy'
    archive = raw / 'amylase_diversity_project-v1.1.zip'
    scratch = ROOT / 'analysis/mixture_null_scratch'
    scratch.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        # Extract to the permitted scratch tree, refusing unsafe archive paths.
        for member in z.infolist():
            p = Path(member.filename)
            if p.is_absolute() or '..' in p.parts:
                raise ValueError('unsafe archive member')
        z.extractall(scratch)
    meta_path = next(scratch.rglob('ancient_metadata_for_supplement.tsv'))
    hap_path = next(scratch.rglob('haplotype_structure_for_supplement.tsv'))
    cn = pd.read_csv(raw/'copy_number_table_for_supplement.tsv', sep='\t')
    meta = pd.read_csv(meta_path, sep='\t')
    hap = pd.read_csv(hap_path, sep='\t')
    if meta['sample'].duplicated().any() or cn['sample'].duplicated().any():
        raise ValueError('ambiguous sample identity')
    ancient = cn[cn.type.eq('Ancient')].merge(meta[['sample','age','region']], on='sample', validate='one_to_one')
    ancient['age'] = pd.to_numeric(ancient.age, errors='coerce')
    ancient['AMY1'] = pd.to_numeric(ancient.AMY1, errors='coerce')
    ancient = ancient.dropna(subset=['age','AMY1','p2']).copy()
    ancient['t'] = -ancient.age/1000
    return ancient, cn, hap


def simulation(t, groups, mu, sigma, effect, reps=1000, seed=SEED):
    """IID conditional CN experiment; group means vary, within-group null flat."""
    g, n, x = design(t, groups)
    sxx = x @ x
    k, N = len(n), len(t)
    if sxx <= 0 or N-k-1 <= 0:
        raise ValueError('no within-group time information')
    rng = np.random.default_rng(seed)
    y = np.asarray(mu)[g, None] + effect*x[:, None] + rng.normal(size=(N,reps))*np.asarray(sigma)[g,None]
    means = np.array([y[g == j].mean(axis=0) for j in range(k)])
    residual = y-means[g]
    betas = x @ residual / sxx
    errors = residual - x[:,None]*betas
    h = 1/n[g] + x*x/sxx
    useful = h < 1-1e-12
    se = np.sqrt(np.sum((x[useful,None]*errors[useful]/(1-h[useful,None]))**2, axis=0))/sxx
    p = 2*stats.t.sf(np.abs(betas/se), N-k-1)
    rejected = int(np.sum(p <= .05))
    ci = stats.binomtest(rejected,reps).proportion_ci()
    # Deliberately uncorrected pooled comparison tests the composition artifact.
    xp = np.asarray(t)-np.mean(t)
    bp = xp @ y/(xp@xp)
    ep = y-y.mean(axis=0)-xp[:,None]*bp
    hp = 1/N+xp*xp/(xp@xp)
    sp = np.sqrt(np.sum((xp[:,None]*ep/(1-hp[:,None]))**2,axis=0))/(xp@xp)
    pooled = float(np.mean(2*stats.t.sf(abs(bp/sp),N-2)<=.05))
    return dict(effect=effect, unit='copies/kyr', seed=seed, replicates=reps,
                rejects=rejected, power=rejected/reps, ci_low=float(ci.low),
                ci_high=float(ci.high), pooled_rejection=pooled)


def amy_analysis(ancient, cn, hap, reps=1000, seed=SEED):
    OUT.mkdir(exist_ok=True)
    groups = []
    for label, d in ancient.groupby('p2'):
        f = fit(d.t, d.AMY1, d.p2)
        groups.append(dict(component=label, mean=float(d.AMY1.mean()),
                           mean_age=float(d.age.mean()), distinct_dates=int(d.age.nunique()),
                           span_kyr=float((d.age.max()-d.age.min())/1000), **f))
    group_df = pd.DataFrame(groups).sort_values('mean_age',ascending=False)
    group_df.to_csv(OUT/'amy_components.tsv',sep='\t',index=False)
    stress = []
    for name, labels in [('pooled',np.repeat('all',len(ancient))),('p2',ancient.p2)]+[
        (f'p2_timebin_{width}kyr',ancient.p2.astype(str)+'_'+np.floor(ancient.age/(width*1000)).astype(int).astype(str))
        for width in [2,1,.5,.25]]+[('singleton',np.arange(len(ancient)).astype(str))]:
        stress.append(dict(model=name, **fit(ancient.t,ancient.AMY1,labels)))
    pd.DataFrame(stress).to_csv(OUT/'amy_dof.tsv',sep='\t',index=False)
    # Independent member means are not fitted to temporal outcomes.
    rng = np.random.default_rng(seed)
    g,n,x = design(ancient.t, ancient.p2)
    means = []; crossmeans = []; sigmas = []
    y = ancient.AMY1.to_numpy(float)
    full = fit(ancient.t, y, ancient.p2)
    for j in range(len(n)):
        indices = np.where(g==j)[0]
        selected = rng.permutation(indices)[:max(1,len(indices)//2)]
        means.append(float(y[indices].mean()))
        crossmeans.append(float(y[selected].mean()))
        sigmas.append(float(np.std(y[indices]-full['beta']*x[indices],ddof=1)) if len(indices)>1 else 0.)
    crossbeta = float(x @ (y-np.asarray(crossmeans)[g])/(x@x))
    curves = [simulation(ancient.t, ancient.p2, means, sigmas, effect,reps,seed) for effect in EFFECTS]
    pd.DataFrame(curves).to_csv(OUT/'amy_power.tsv',sep='\t',index=False)
    refined=[simulation(ancient.t,ancient.p2,means,sigmas,effect,reps,seed) for effect in [.35,.4,.45]]
    pd.DataFrame(refined).to_csv(OUT/'amy_power_refined.tsv',sep='\t',index=False)
    amy_decomposition(ancient,seed)
    # Slot-weighted carrier diploid CN, explicitly not haploid copy-number means.
    merged = hap.merge(cn[['sample','AMY1']],on='sample',validate='one_to_one')
    slots = pd.concat([merged.rename(columns={col:'haplotype'}).drop(columns=[other])
                       for col,other in [('H1','H2'),('H2','H1')]],ignore_index=True)
    slots.groupby('haplotype').agg(slot_n=('AMY1','size'),carrier_diploid_mean=('AMY1','mean')).to_csv(OUT/'haplotype_means.tsv',sep='\t')
    sets = {label:sorted(set(slots.loc[slots.type.eq(label),'haplotype'].dropna())) for label in ['Ancient','Modern']}
    modern = slots[slots.type.eq('Modern')].copy()
    modern['H5'] = modern.haplotype.eq('H5')
    for col in ['p1','p2']:
        modern.groupby(col).agg(slots=('H5','size'),H5_frequency=('H5','mean')).to_csv(OUT/f'h5_by_{col}.tsv',sep='\t')
    audit = dict(n_ancient=len(ancient),regions=ancient.region.value_counts().to_dict(),
                 ancient_haplotypes=sets['Ancient'],modern_haplotypes=sets['Modern'],
                 ancient_only=sorted(set(sets['Ancient'])-set(sets['Modern'])),
                 full_beta=full['beta'],crossfit_baseline_beta=crossbeta,
                 crossfit_baseline_max_slope_difference=abs(crossbeta-full['beta']),
                 grid_mde80=next((c['effect'] for c in curves if c['power']>=.8),None),
                 refined_grid_mde80=next((c['effect'] for c in sorted(curves+refined,key=lambda c:c['effect']) if c['power']>=.8),None),
                 conditional_null_gate=curves[0]['ci_low']<=.05<=curves[0]['ci_high'],
                 seed=seed,replicates=reps)
    (OUT/'amy_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes = plt.subplots(1,3,figsize=(15,4.5))
    for label,d in ancient.groupby('p2'):
        axes[0].scatter(d.t,d.AMY1,s=9,alpha=.35,label=label)
    axes[0].set(xlabel='Forward time (-BP/1000)',ylabel='AMY1 copies',title='Measured labeled individuals')
    axes[0].legend(fontsize=6)
    axes[1].errorbar(group_df.beta,range(len(group_df)),xerr=1.96*group_df.se,fmt='o')
    axes[1].set_yticks(range(len(group_df)),group_df.component,fontsize=7)
    axes[1].axvline(0,color='gray'); axes[1].set(xlabel='Within-label copies/kyr',title='Slopes ±1.96 HC3 SE')
    axes[2].plot(EFFECTS,[c['power'] for c in curves],'o-',label='Within-group residual')
    axes[2].plot(EFFECTS,[c['pooled_rejection'] for c in curves],'--',label='Uncorrected pooled')
    axes[2].axhline(.8,color='gray',ls=':'); axes[2].set(xlabel='Injected copies/kyr',ylabel='Rejection fraction',ylim=(0,1.03),title='Conditional power and specificity')
    axes[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(ROOT/'fig/amy_mixture_null.png',dpi=170); plt.close(fig)
    return audit


def amy_decomposition(ancient, seed=SEED):
    """Descriptive composition slopes; baseline holdout is not causal attribution."""
    t=ancient.t.to_numpy(); y=ancient.AMY1.to_numpy(float)
    g,n,x=design(t,ancient.p2); mu=np.bincount(g,weights=y)/n
    xp=t-t.mean(); total=float(xp@y/(xp@xp)); comp=float(xp@mu[g]/(xp@xp))
    full=fit(t,y,ancient.p2); rng=np.random.default_rng(seed)
    cross=[]; held=[]
    for j in range(len(n)):
        ix=np.where(g==j)[0]; perm=rng.permutation(ix); split=max(1,len(ix)//2)
        cross.append(float(y[perm[:split]].mean())); held.extend(perm[split:])
    cross=np.array(cross); held=np.array(held,int)
    result=dict(pooled_beta=total,composition_beta=comp,composition_slope_fraction=comp/total,
                within_contribution=full['beta']*(x@x)/(xp@xp),within_time_information_fraction=(x@x)/(xp@xp),
                crossfit_composition_beta=float(xp@cross[g]/(xp@xp)),
                heldout_baseline_RMSE=float(np.sqrt(np.mean((y[held]-cross[g[held]])**2))),heldout_baseline_n=len(held))
    (OUT/'amy_decomposition.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def dof_power(ancient, reps=1000, seed=SEED):
    """Same injected dataset, increasingly flexible pre-fixed nuisance labels."""
    t=ancient.t.to_numpy(); y=ancient.AMY1.to_numpy(float)
    g,n,x=design(t,ancient.p2)
    mu=np.bincount(g,weights=y)/n
    f=fit(t,y,ancient.p2)
    sigma=np.array([np.std(y[g==j]-f['beta']*x[g==j],ddof=1) if n[j]>1 else 0. for j in range(len(n))])
    rng=np.random.default_rng(seed)
    noise=rng.normal(size=(len(y),reps))*sigma[g,None]
    partitions=[('pooled',np.repeat('all',len(y))),('p2',ancient.p2)]+[
        (f'p2_timebin_{w}kyr',ancient.p2.astype(str)+'_'+np.floor(ancient.age/(w*1000)).astype(int).astype(str))
        for w in [2,1,.5,.25]]+[('singleton',np.arange(len(y)).astype(str))]
    rows=[]
    for effect in [0,.5,2]:
        Y=mu[g,None]+effect*x[:,None]+noise
        for label,groups in partitions:
            gg,nn,xx=design(t,groups); N=len(t); K=len(nn); sxx=xx@xx; df=N-K-1
            if sxx<1e-12 or df<=0:
                rows.append(dict(effect=effect,model=label,k=K,df=df,identified=False,power=None))
                continue
            means=np.array([Y[gg==j].mean(axis=0) for j in range(K)])
            b=xx@(Y-means[gg])/sxx; e=Y-means[gg]-xx[:,None]*b
            h=1/nn[gg]+xx*xx/sxx; ok=h<1-1e-12
            se=np.sqrt(np.sum((xx[ok,None]*e[ok]/(1-h[ok,None]))**2,axis=0))/sxx
            p=2*stats.t.sf(abs(b/se),df)
            count=int((p<=.05).sum()); ci=stats.binomtest(count,reps).proportion_ci()
            rows.append(dict(effect=effect,model=label,k=K,df=df,identified=True,power=count/reps,
                             ci_low=float(ci.low),ci_high=float(ci.high),seed=seed,replicates=reps))
    pd.DataFrame(rows).to_csv(OUT/'amy_dof_power.tsv',sep='\t',index=False)
    return rows


def strict_amy_feasibility(ancient):
    """Independent exact-token join; no calls to previous result-writing helper."""
    import csv
    raw=ROOT/'data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.anno'
    tokens={}
    with raw.open() as fh:
        reader=csv.reader(fh,delimiter='\t'); next(reader)
        for row in reader:
            for tok in {row[0].split('.')[0],row[2],row[3]}-{'','..','.'}:
                tokens.setdefault(tok,set()).add(row[2])
    members=pd.read_csv(ROOT/'analysis/deme_rebuild/individuals.tsv',sep='\t')
    index=members.set_index('biological_id')
    joined=[]; ambiguous=0; unmatched=0
    for row in ancient.itertuples():
        ids=tokens.get(row.sample,set())
        if len(ids)>1:
            ambiguous+=1; continue
        if not ids or next(iter(ids)) not in index.index:
            unmatched+=1; continue
        a=index.loc[next(iter(ids))]
        joined.append(dict(sample=row.sample,AMY1=row.AMY1,p2=row.p2,archive_age=row.age,
                           biological_id=a.name,deme_id=int(a.deme_id),node_id=int(a.node_id),
                           date_BP=float(a.date_BP)))
    j=pd.DataFrame(joined)
    if j.biological_id.duplicated().any():
        raise ValueError('multiple AMY observations of one biological identity')
    j.to_csv(OUT/'amy_strict_joined.tsv',sep='\t',index=False)
    rows=[]; node_rows=[]
    for d,part in j[j.deme_id.ge(0)].groupby('deme_id'):
        nodes=part.groupby('node_id').agg(n=('AMY1','size'),mean_CN=('AMY1','mean'),date_BP=('date_BP','mean')).reset_index()
        nodes['deme_id']=d; node_rows.append(nodes)
        supported=nodes[nodes.n.ge(5)]
        eligible=part[part.node_id.isin(supported.node_id)]
        span=float(np.ptp(part.date_BP)/1000)
        supported_span=float(np.ptp(eligible.date_BP)/1000) if len(eligible) else 0.
        usable=len(supported)>=3 and supported_span>=2
        sf=fit(-eligible.date_BP/1000,eligible.AMY1,np.repeat('all',len(eligible))) if usable else {}
        rows.append(dict(deme_id=int(d),n=len(part),nodes=len(nodes),span_kyr=span,
                         supported_nodes=len(supported),supported_n=len(eligible),supported_span_kyr=supported_span,
                         strict_usable=usable,**{f'supported_fit_{k}':v for k,v in sf.items()}))
    d=pd.DataFrame(rows); d.to_csv(OUT/'amy_strict_demes.tsv',sep='\t',index=False)
    if node_rows: pd.concat(node_rows).to_csv(OUT/'amy_strict_nodes.tsv',sep='\t',index=False)
    usable=d[d.strict_usable]
    f=fit(-j[j.deme_id.ge(0)].date_BP/1000,j[j.deme_id.ge(0)].AMY1,j[j.deme_id.ge(0)].deme_id.astype(str))
    result=dict(release_ancient_n=len(ancient),linked_rebuilt_members=len(j),ambiguous_tokens=ambiguous,
                unmatched_or_not_retained=unmatched,carriers_in_qualifying_paths=int(j.deme_id.ge(0).sum()),
                paths_with_carriers=len(d),strict_usable_paths=len(usable),strict_replication_possible=len(usable)>=2,
                exploratory_all_carrier_common_fit=f)
    (OUT/'amy_strict_feasibility.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--simulate',action='store_true')
    ap.add_argument('--effect',type=float,default=.5)
    ap.add_argument('--power-curve',action='store_true')
    ap.add_argument('--dof-power',action='store_true')
    ap.add_argument('--strict-feasibility',action='store_true')
    ap.add_argument('--reps',type=int,default=1000)
    ap.add_argument('--seed',type=int,default=SEED)
    args=ap.parse_args()
    OUT.mkdir(exist_ok=True)
    ancient,cn,hap=load_amy()
    if args.dof_power:
        result=dof_power(ancient,args.reps,args.seed)
    elif args.strict_feasibility:
        result=strict_amy_feasibility(ancient)
    elif args.simulate or args.power_curve:
        labels,g = np.unique(ancient.p2.astype(str),return_inverse=True)
        means=[ancient.AMY1.to_numpy()[g==j].mean() for j in range(len(labels))]
        f=fit(ancient.t,ancient.AMY1,ancient.p2)
        _,_,x=design(ancient.t,ancient.p2)
        sigma=[np.std(ancient.AMY1.to_numpy()[g==j]-f['beta']*x[g==j],ddof=1) if np.sum(g==j)>1 else 0. for j in range(len(labels))]
        effects=EFFECTS if args.power_curve else [args.effect]
        result=[simulation(ancient.t,ancient.p2,means,sigma,e,args.reps,args.seed) for e in effects]
    else:
        result=amy_analysis(ancient,cn,hap,args.reps,args.seed)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
