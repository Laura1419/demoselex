#!/usr/bin/env python3
"""Spatially anchored candidate deme paths; prereg: docs/deme-rebuild-prereg.md."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from analysis.load_aadr import BUNDLE, MAIN_PREFIX, ANNO_NAME, load_bundle, iter_snp
from analysis.traj_relatedness import relatedness_from_accum
from scripts import provenance

OUT = Path('analysis/deme_rebuild')
SEED = 20261008


def recorder(slug):
    # The user explicitly forbids ANY git invocation, including provenance probes.
    provenance.get_git_info = lambda: {'commit': 'not-probed:no-git-rule',
                                      'branch': 'not-probed', 'dirty': 'unknown'}
    return provenance.RunRecorder(slug)


def build_paths(meta, pcs):
    """No component merging: one-to-one, adjacent-bin, anchor-bounded chains."""
    nodes, paths = [], []
    for cell, part in meta.groupby(['lat_cell', 'lon_cell'], sort=True):
        local = []
        for time_bin, chunk in part.groupby('time_bin', sort=True):
            ix = chunk.index.to_numpy()
            if len(ix) < 8:
                continue
            k = max(1, min(4, len(ix)//20))
            labels = (KMeans(n_clusters=k, n_init=5, random_state=SEED).fit_predict(pcs[ix])
                      if k > 1 else np.zeros(len(ix), dtype=int))
            for label in sorted(set(labels)):
                members = ix[labels == label]
                if len(members) < 8:
                    continue
                centroid = pcs[members].mean(axis=0)
                radius = float(np.sqrt(np.mean(np.sum((pcs[members]-centroid)**2, axis=1))))
                node = dict(node_id=len(nodes), cell=cell, time_bin=int(time_bin),
                            members=members, centroid=centroid, radius=radius,
                            date_BP=float(meta.loc[members, 'date_BP'].mean()))
                nodes.append(node)
                local.append(node)
        if not local:
            continue
        radius = float(np.median([n['radius'] for n in local]))
        previous, active = [], {}
        for time_bin in sorted({n['time_bin'] for n in local}):
            current = [n for n in local if n['time_bin'] == time_bin]
            next_active = {}
            if previous and previous[0]['time_bin'] == time_bin-1:
                distance = np.array([[np.linalg.norm(a['centroid']-b['centroid'])
                                     for b in current] for a in previous])
                allowed = np.zeros_like(distance, dtype=bool)
                for i, a in enumerate(previous):
                    path = active[a['node_id']]
                    for j, b in enumerate(current):
                        dates = meta.loc[np.concatenate([n['members'] for n in path]+[b['members']]), 'date_BP']
                        allowed[i,j] = (distance[i,j] <= 1.5*radius and
                            np.linalg.norm(path[0]['centroid']-b['centroid']) <= 2*radius and
                            dates.max()-dates.min() <= 8000)
                rr, cc = linear_sum_assignment(np.where(allowed, distance, 1e12))
                for i, j in zip(rr, cc):
                    if allowed[i,j]:
                        path = active[previous[i]['node_id']]
                        path.append(current[j])
                        next_active[current[j]['node_id']] = path
            for b in current:
                if b['node_id'] not in next_active:
                    path = [b]
                    paths.append(path)
                    next_active[b['node_id']] = path
            previous, active = current, next_active
    return nodes, paths


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rec = recorder('deme-rebuild')
    rec.set_seed(SEED)
    for f in ['analysis/deme_rebuild.py', 'analysis/load_aadr.py',
              'docs/deme-rebuild-prereg.md', f'{BUNDLE}/{ANNO_NAME}',
              f'{BUNDLE}/{MAIN_PREFIX}.ind', f'{BUNDLE}/{MAIN_PREFIX}.snp']:
        rec.record_input(f)
    gm, individuals, anno = load_bundle()
    rec.record_parameter('raw_genotype', {'path': gm.path, 'bytes': Path(gm.path).stat().st_size})
    records = []
    for row, ind in enumerate(individuals):
        a = anno.get(ind.id)
        if (a is None or a.is_modern or not a.has_date or not a.has_coords or
            not (0 < a.date_bp <= 14000) or not np.isfinite([a.date_bp,a.latitude,a.longitude]).all() or
            any(x in a.assessment.upper() for x in ['FAIL', 'CRITICAL'])):
            continue
        records.append(dict(individual=ind.id, biological_id=a.individual_id,
                            row=row, date_BP=a.date_bp, lat=a.latitude, lon=a.longitude,
                            published_group=a.group_id))
    meta = pd.DataFrame(records)
    candidate = []
    for i, s in enumerate(iter_snp(f'{BUNDLE}/{MAIN_PREFIX}.snp')):
        if s[1] in {str(c) for c in range(2,23) if c != 6}:
            candidate.append(i)
    rng = np.random.default_rng(SEED)
    chosen = rng.choice(candidate, size=7800, replace=False)
    embedding_sites, control_sites = np.sort(chosen[:6000]), chosen[6000:]
    codes = np.stack([gm.individual_genotypes(i, embedding_sites) for i in meta.row])
    meta['presence'] = (codes != 3).mean(axis=1)
    before = len(meta)
    keep = meta.presence >= .2
    dedup = (meta.loc[keep].sort_values(['presence','row'], ascending=[False,True])
             .drop_duplicates('biological_id').sort_values('row'))
    codes = codes[dedup.index]
    meta = dedup.reset_index(drop=True)
    X = np.where(codes == 3, np.nan, codes/2).astype(np.float32)
    del codes
    freq = np.nanmean(X, axis=0)
    good = np.isfinite(freq) & (freq > .01) & (freq < .99)
    X = X[:,good]
    freq = freq[good]
    X = np.where(np.isnan(X), freq, X)
    X = (X-freq)/np.sqrt(freq*(1-freq))
    pca = PCA(n_components=8, svd_solver='randomized', random_state=SEED)
    pcs = pca.fit_transform(X)
    del X
    meta['lat_cell'] = np.floor(meta.lat/5).astype(int)
    meta['lon_cell'] = np.floor(meta.lon/5).astype(int)
    meta['time_bin'] = np.floor(meta.date_BP/1000).astype(int)
    nodes, paths = build_paths(meta, pcs)
    summary, node_table = [], []
    meta['deme_id'] = -1
    meta['node_id'] = -1
    qualifying = []
    for path_id, path in enumerate(paths):
        members = np.concatenate([n['members'] for n in path])
        dates = meta.loc[members,'date_BP']
        span = float(dates.max()-dates.min())
        qualifies = len(path)>=3 and len(members)>=30 and span>=2000
        row = dict(path_id=path_id, qualifies=qualifies, n=len(members), n_nodes=len(path),
                   min_node_n=min(len(n['members']) for n in path),
                   date_min=float(dates.min()), date_max=float(dates.max()), span=span,
                   lat_cell=path[0]['cell'][0],lon_cell=path[0]['cell'][1])
        summary.append(row)
        if not qualifies:
            continue
        deme_id = len(qualifying)
        row['deme_id'] = deme_id
        qualifying.append(path)
        for n in path:
            meta.loc[n['members'],'deme_id'] = deme_id
            meta.loc[n['members'],'node_id'] = n['node_id']
            node_table.append(dict(deme_id=deme_id,node_id=n['node_id'],
                time_bin=n['time_bin'], date_BP=n['date_BP'], n=len(n['members']),
                **{f'PC{i+1}': float(v) for i,v in enumerate(n['centroid'])}))
    paths_df = pd.DataFrame(summary)
    qualified_df = paths_df.loc[paths_df.qualifies].copy()
    if len(qualified_df):
        qualified_df['deme_id'] = qualified_df.deme_id.astype(int)
    meta.to_csv(OUT/'individuals.tsv',sep='\t',index=False)
    paths_df.to_csv(OUT/'paths.tsv',sep='\t',index=False)
    pd.DataFrame(node_table).to_csv(OUT/'nodes.tsv',sep='\t',index=False)
    np.savez_compressed(OUT/'panels.npz', embedding_sites=embedding_sites, control_sites=control_sites)
    if qualifying:
        sums, counts = [], []
        for path in qualifying:
            members = np.concatenate([n['members'] for n in path])
            g = np.stack([gm.individual_genotypes(int(meta.loc[i,'row']), embedding_sites)
                          for i in members])
            sums.append(np.where(g!=3,g,0).sum(axis=0))
            counts.append((g!=3).sum(axis=0))
        mat, used = relatedness_from_accum(np.arange(len(qualifying)),np.array(sums,dtype=float),np.array(counts))
        mat.to_csv(OUT/'ancestry_overlap.tsv',sep='\t')
        neff = float(len(mat)**2/mat.to_numpy().sum())
        off = mat.to_numpy()[np.triu_indices(len(mat),1)]
    else:
        used, neff, off = 0, None, np.array([])
    metrics = dict(timestamp=pd.Timestamp.now(tz='UTC').isoformat(), eligible_before_coverage=before,
        after_coverage=int(keep.sum()), after_dedup=len(meta), nodes=len(nodes), paths=len(paths),
        qualifying_demes=len(qualifying), assigned=int((meta.deme_id>=0).sum()),
        ancestry_overlap_sites=used, ancestry_neff_diagnostic=neff,
        max_ancestry_overlap=float(off.max()) if len(off) else None,
        median_ancestry_overlap=float(np.median(off)) if len(off) else None,
        PCA_variance=pca.explained_variance_ratio_.tolist(),
        control_locus_pass_fraction=None, calibration='not yet run; see directional_consistency.py')
    (OUT/'metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    for f in OUT.iterdir():
        if f.is_file(): rec.record_output(f)
    rec.save()
    print(json.dumps(metrics,indent=2),flush=True)
    print(qualified_df.to_string(index=False),flush=True)


if __name__ == '__main__':
    main()
