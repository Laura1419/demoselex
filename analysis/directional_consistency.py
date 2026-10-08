#!/usr/bin/env python3
"""Held-out drift-plus-jump calibration and AMY feasibility (preregistered)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.stats import beta, norm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from analysis.deme_rebuild import OUT, SEED, recorder
from analysis.load_aadr import BUNDLE, MAIN_PREFIX, GenotypeMatrix
from analysis.amy import load_copy_number_table, _aadr_rows, ANCIENT_TYPES


def slope_weights(time, count):
    time, count = np.asarray(time,float), np.asarray(count,float)
    centered = time - np.average(time,weights=count)
    return count*centered / np.sum(count*centered**2)


def process_variance(weights, sampling, increment_variance):
    """Cumulative innovations induce full path covariance; sampling is diagonal."""
    innovation_weights = np.cumsum(weights[::-1])[::-1][1:]
    return np.sum(weights**2*sampling)+np.sum(innovation_weights**2*increment_variance)


def joint_score(z, correlation):
    return np.sum(z,axis=-1)/np.sqrt(np.sum(correlation))


def conjunction_score(z):
    """Weakest replicate in either direction, not an outlier-dominated mean."""
    z = np.asarray(z,float)
    return np.maximum(np.min(z,axis=-1),np.min(-z,axis=-1))


def interval(k,n):
    if not n: return [None,None]
    return [float(beta.ppf(.025,k,n-k+1)) if k else 0.,
            float(beta.ppf(.975,k+1,n-k)) if k<n else 1.]


def amy_feasibility(meta):
    # Use existing CN loader/token definitions, but reject ambiguous biological IDs.
    token_ids = {}
    raw = _aadr_rows()
    for row in raw:
        tokens = {row[0].split('.')[0]} | {row[c].strip() for c in (2,3)
                                          if row[c].strip() not in ('','..','.')}
        for token in tokens:
            token_ids.setdefault(token,set()).add(row[2])
    by_id = meta.set_index('biological_id')
    cn = [r for r in load_copy_number_table() if r['type'] in ANCIENT_TYPES]
    joined = []
    ambiguous = numeric = token_matched = 0
    for r in cn:
        ids = token_ids.get(r['sample'],set())
        if len(ids)>1:
            ambiguous += 1
            continue
        if not ids: continue
        token_matched += 1
        biological_id = next(iter(ids))
        # Validate real schema columns, not truthiness of unrelated metadata strings.
        matches = [a for a in raw if a[2]==biological_id]
        valid = False
        for a in matches:
            try: valid |= bool(np.isfinite([float(a[c]) for c in (10,17,18)]).all())
            except ValueError: pass
        numeric += int(valid)
        if biological_id not in by_id.index: continue
        a = by_id.loc[biological_id]
        joined.append(dict(sample=r['sample'],AMY1=int(r['AMY1']),source=r['source'],
                           biological_id=biological_id,individual=a.individual,
                           deme_id=int(a.deme_id),node_id=int(a.node_id),date_BP=float(a.date_BP)))
    j = pd.DataFrame(joined)
    j.to_csv(OUT/'amy_joined.tsv',sep='\t',index=False)
    rows = []
    if len(j):
        for deme_id, d in j.loc[j.deme_id>=0].groupby('deme_id',sort=True):
            node = d.groupby('node_id').agg(n=('AMY1','size'),mean_CN=('AMY1','mean'),
                                          date_BP=('date_BP','mean')).reset_index()
            node['deme_id'] = deme_id
            rows.append((d,node))
    summary, all_nodes = [], []
    for d, node in rows:
        all_nodes.append(node)
        supported = node.loc[node.n>=5]
        support_samples = d.loc[d.node_id.isin(supported.node_id)]
        span = float(d.date_BP.max()-d.date_BP.min())
        supported_span = (float(support_samples.date_BP.max()-support_samples.date_BP.min())
                          if len(support_samples) else 0.)
        usable = len(supported)>=3 and supported_span>=2000
        descriptive_slope = None
        if len(node)>=3 and span>=2000:
            w = slope_weights(-node.date_BP.to_numpy()/1000,node.n.to_numpy())
            descriptive_slope = float(w @ node.mean_CN.to_numpy())
        supported_slope = None
        oldest_mean = newest_mean = None
        if usable:
            supported = supported.sort_values('date_BP',ascending=False)
            w = slope_weights(-supported.date_BP.to_numpy()/1000,supported.n.to_numpy())
            supported_slope = float(w @ supported.mean_CN.to_numpy())
            oldest_mean = float(supported.mean_CN.iloc[0])
            newest_mean = float(supported.mean_CN.iloc[-1])
        summary.append(dict(deme_id=int(d.deme_id.iloc[0]),n=len(d),n_nodes=len(node),
                            min_node_n=int(node.n.min()),date_min=float(d.date_BP.min()),
                            date_max=float(d.date_BP.max()),span=span,n_supported_nodes=len(supported),
                            supported_span=supported_span,qualifies_CN=usable,
                            descriptive_slope_copies_per_kyr=descriptive_slope,
                            supported_slope_copies_per_kyr=supported_slope,
                            oldest_supported_mean_CN=oldest_mean,newest_supported_mean_CN=newest_mean))
    pd.DataFrame(summary).to_csv(OUT/'amy_demes.tsv',sep='\t',index=False)
    if all_nodes: pd.concat(all_nodes).to_csv(OUT/'amy_nodes.tsv',sep='\t',index=False)
    qualifying = [r for r in summary if r['qualifies_CN']]
    eligible = [r for r in summary if r['n']>=5]
    cells = {}
    for r in eligible:
        cell = meta.loc[meta.deme_id==r['deme_id'],['lat_cell','lon_cell']].iloc[0]
        cells.setdefault(f'{int(cell.lat_cell)},{int(cell.lon_cell)}',[]).append(r)
    successions = []
    for cell, series in cells.items():
        series = sorted(series,key=lambda r: -r['date_max'])
        if len(series)>=3 and all(a['date_min']>b['date_max'] for a,b in zip(series,series[1:])):
            successions.append(dict(cell=cell,demes=[r['deme_id'] for r in series]))
    return dict(n_ancient=len(cn),token_matched=token_matched,ambiguous=ambiguous,
                numeric_date_coordinates=numeric,retained_rebuild_samples=len(j),
                in_qualifying_demes=int((j.deme_id>=0).sum()) if len(j) else 0,
                demes_with_carriers=len(summary),demes_with_three_nodes_and_2000yr=
                sum(r['n_nodes']>=3 and r['span']>=2000 for r in summary),
                qualifying_CN_demes=len(qualifying),demes=summary,
                ordered_deme_paths_with_5_carriers=len(eligible),
                max_CN_paths_per_cell=max((len(v) for v in cells.values()),default=0),
                qualifying_nonoverlapping_successions=successions,
                drift_jump_ordering_probability=None,
                CN_calibration='not identifiable: no independently calibrated CN process or CN controls')


def main():
    rec = recorder('directional-consistency')
    rec.set_seed(SEED)
    for f in ['analysis/directional_consistency.py','docs/deme-rebuild-prereg.md',
              'analysis/amy.py','data/raw/amy/copy_number_table_for_supplement.tsv',
              str(OUT/'individuals.tsv'),str(OUT/'nodes.tsv'),str(OUT/'paths.tsv'),str(OUT/'panels.npz')]:
        rec.record_input(f)
    meta = pd.read_csv(OUT/'individuals.tsv',sep='\t')
    paths = pd.read_csv(OUT/'paths.tsv',sep='\t')
    panel = (paths.loc[paths.qualifies & (paths.min_node_n>=20)]
             .sort_values(['min_node_n','span','path_id'],ascending=[False,False,True]).head(6))
    ids = panel.deme_id.astype(int).tolist()
    metrics = dict(timestamp=pd.Timestamp.now(tz='UTC').isoformat(),process_panel_demes=ids,
                   calibration='infeasible',control_locus_pass_fraction=None)
    print('Prespecified process panel:',ids,flush=True)
    if len(ids)>=2:
        nodes = pd.read_csv(OUT/'nodes.tsv',sep='\t')
        sites = np.load(OUT/'panels.npz')['control_sites']
        gm = GenotypeMatrix(f'{BUNDLE}/{MAIN_PREFIX}.geno')
        frequencies, counts, times, distances = [], [], [], []
        for deme_id in ids:
            ns = nodes.loc[nodes.deme_id==deme_id].sort_values('date_BP',ascending=False)
            fs, cs = [], []
            for node_id in ns.node_id:
                rr = meta.loc[meta.node_id==node_id,'row']
                g = np.stack([gm.individual_genotypes(int(r),sites) for r in rr])
                count = (g!=3).sum(axis=0)
                fs.append(np.divide(np.where(g!=3,g,0).sum(axis=0),2*count,
                                    out=np.zeros(len(sites),dtype=float),where=count>0))
                cs.append(count)
            frequencies.append(np.array(fs).T)
            counts.append(np.array(cs).T)
            times.append(-ns.date_BP.to_numpy()/1000)
            distances.append(np.sum(np.diff(ns[[f'PC{i}' for i in range(1,9)]].to_numpy(),axis=0)**2,axis=1))
        complete = np.ones(len(sites),dtype=bool)
        for c in counts: complete &= (c>=8).all(axis=1)
        af = np.mean([f.mean(axis=1) for f in frequencies],axis=0)
        complete &= (af>=.1)&(af<=.9)
        train = np.flatnonzero(complete & (np.arange(len(sites))<900))
        validation = np.flatnonzero(complete & (np.arange(len(sites))>=900))
        print('Complete controls train/validation:',len(train),len(validation),flush=True)
        metrics.update(training_controls=len(train),validation_controls=len(validation))
        if len(train)>=100 and len(validation)>=100:
            design, targets = [], []
            for f,c,t,d2 in zip(frequencies,counts,times,distances):
                base = np.clip(f[train].mean(axis=1),.01,.99)
                scale = base*(1-base)
                diff = np.diff(f[train],axis=1)
                target = diff**2/scale[:,None] - (1/c[train,:-1]+1/c[train,1:])
                features = np.column_stack([np.diff(t)*1000/29,d2,np.ones(len(d2))])
                design.append(np.tile(features,(len(train),1)))
                targets.append(target.ravel())
            X, y = np.vstack(design),np.concatenate(targets)
            # Column scaling changes conditioning, not the NNLS model.
            scaleX = np.sqrt(np.mean(X**2,axis=0))
            scaleX[scaleX==0] = 1
            coefficients, _ = nnls(X/scaleX,y)
            coefficients /= scaleX
            # Equivalent sufficient-statistic regression, resampling whole loci jointly.
            unique_design = np.vstack([np.column_stack([np.diff(t)*1000/29,d2,np.ones(len(d2))])
                                      for t,d2 in zip(times,distances)])
            locus_targets = np.hstack([v.reshape(len(train),-1) for v in targets])
            rng = np.random.default_rng(SEED+1)
            boot = []
            for _ in range(500):
                pick = rng.integers(0,len(train),len(train))
                fitted, _ = nnls(unique_design/scaleX,locus_targets[pick].mean(axis=0))
                boot.append(fitted/scaleX)
            boot = np.array(boot)
            aq = np.quantile(boot[:,0],[.025,.5,.975])
            metrics['drift_fit_bootstrap'] = dict(n=500,zero_boundary_fraction=float((boot[:,0]==0).mean()),
                drift_quantiles_025_50_975=aq.tolist(),
                Ne_quantiles_025_50_975=[float(1/(2*a)) if a>0 else 'infinity' for a in aq[::-1]])
            z, slopes, changes = [], [], []
            for f,c,t,d2 in zip(frequencies,counts,times,distances):
                q = np.column_stack([np.diff(t)*1000/29,d2,np.ones(len(d2))]) @ coefficients
                zz, ss, dd = [], [], []
                for index in range(len(sites)):
                    if not complete[index]:
                        zz.append(np.nan); ss.append(np.nan); dd.append(np.nan)
                        continue
                    weights = slope_weights(t,c[index])
                    base = np.clip(f[index].mean(),.01,.99)
                    variance = process_variance(weights,base*(1-base)/np.maximum(c[index],1),base*(1-base)*q)
                    slope = float(weights @ f[index])
                    zz.append(slope/np.sqrt(variance))
                    ss.append(slope)
                    dd.append(slope*(t.max()-t.min()))
                z.append(zz); slopes.append(ss); changes.append(dd)
            z = np.array(z).T
            correlation = .5*np.corrcoef(z[train],rowvar=False)+.5*np.eye(len(ids))
            pd.DataFrame(correlation,index=ids,columns=ids).to_csv(OUT/'slope_correlation.tsv',sep='\t')
            score = joint_score(z,correlation)
            p = 2*norm.sf(abs(score))
            null_z = np.random.default_rng(SEED+2).multivariate_normal(
                np.zeros(len(ids)),correlation,size=100000)
            null_conjunction = np.sort(conjunction_score(null_z))
            conjunction = conjunction_score(z)
            p_conjunction = (1+len(null_conjunction)-np.searchsorted(
                null_conjunction,conjunction,side='left'))/(1+len(null_conjunction))
            p_conjunction = np.where(complete,p_conjunction,np.nan)
            threshold = float(np.quantile(null_conjunction,.95))
            conjunction_k = int((p_conjunction[validation]<=.05).sum())
            conjunction_ci = interval(conjunction_k,len(validation))
            metrics['strict_directional_conjunction'] = dict(null_draws=100000,
                control_pass_count=conjunction_k,control_pass_fraction=conjunction_k/len(validation),
                control_pass_95CI=conjunction_ci,
                calibration='passed' if conjunction_ci[0]<=.05<=conjunction_ci[1] else 'FAILED: interpretation blocked',
                null_95pct_threshold=threshold,
                conditional_power={str(effect):float((conjunction_score(null_z+effect)>threshold).mean())
                                   for effect in [.25,.5,1.,2.]})
            df = pd.DataFrame(dict(site_index=sites,AF=af,complete=complete,
                                   split=np.where(np.arange(len(sites))<900,'train','validation'),
                                   score=score,p_process=p,n_positive=(z>0).sum(axis=1),
                                   conjunction=conjunction,p_conjunction=p_conjunction))
            for i,deme_id in enumerate(ids):
                df[f'slope_{deme_id}']=np.array(slopes[i])
                df[f'change_{deme_id}']=np.array(changes[i])
                df[f'z_{deme_id}']=z[:,i]
            df.to_csv(OUT/'control_loci.tsv',sep='\t',index=False)
            n = len(validation); k = int((p[validation]<=.05).sum()); ci = interval(k,n)
            effective = len(ids)**2/float(correlation.sum())
            power = {str(effect): float(norm.sf(1.959963984540054-effect*np.sqrt(effective))+
                                       norm.cdf(-1.959963984540054-effect*np.sqrt(effective)))
                     for effect in [.25,.5,1.,2.]}
            metrics.update(coefficients=dict(drift_per_generation=float(coefficients[0]),
                jump_per_PC_distance_squared=float(coefficients[1]),jump_per_transition=float(coefficients[2])),
                Ne_estimate=float(1/(2*coefficients[0])) if coefficients[0]>0 else None,
                scaled_design_condition=float(np.linalg.cond(X/scaleX)),
                slope_information_equivalent_D=effective,conditional_power=power,
                control_locus_pass_count=k,control_locus_pass_fraction=k/n,control_pass_95CI=ci,
                calibration='passed' if ci[0]<=.05<=ci[1] else 'FAILED: biological interpretation blocked',
                validation_mean_score=float(score[validation].mean()),
                validation_sd_score=float(score[validation].std(ddof=1)))
    # Feasibility/descriptive CN only; no interpretation bypasses SNP or CN gates.
    metrics['AMY'] = amy_feasibility(meta)
    (OUT/'directional_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    for name in ['directional_metrics.json','control_loci.tsv','slope_correlation.tsv',
                 'amy_joined.tsv','amy_demes.tsv','amy_nodes.tsv']:
        f = OUT/name
        if f.exists(): rec.record_output(f)
    rec.save()
    print(json.dumps(metrics,indent=2),flush=True)


if __name__ == '__main__':
    main()
