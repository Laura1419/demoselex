#!/usr/bin/env python3
"""Full-panel unphased window backgrounds; docs/haplotype-background-prereg.md."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import re
import resource
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, t as student_t
from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.load_akbari_geno import AkbariGeno, paths, parse_psam

OUT = Path('data/derived/haplotype_background')
PRE = Path('docs/haplotype-background-prereg.md')
SEED = 20261009
DIM = 48
KS = (2, 4, 6, 8)
LOCI = {'AMY': (1, 103915000, 104255000),
        'LCT': (2, 136608646-25000, 136608646+25000),
        'FADS': (11, 61569830-25000, 61569830+25000),
        'SLC45A2': (5, 33951693-25000, 33951693+25000),
        'SLC24A5': (15, 48426484-25000, 48426484+25000)}
ANNO = Path('data/raw/aadr-v66.p1/v66.p1_1240K.aadr.PUB.anno')
ARCHIVE = Path('data/raw/amy/amylase_diversity_project-v1.1.zip')
ANNOTATION = Path('data/raw/annotation/hg19_refGene.txt.gz')
FOCAL_GENES = {'AMY': ['AMY1A', 'AMY1B', 'AMY1C', 'AMY2A', 'AMY2B'],
               'LCT': ['LCT'], 'FADS': ['FADS1', 'FADS2', 'FADS3'],
               'SLC45A2': ['SLC45A2'], 'SLC24A5': ['SLC24A5']}
EXCLUSIONS = None
EXCLUSION_AUDIT = None


def initialize_intervals():
    """Annotation-derived, transcript-complete cores with the unchanged 25kb buffer."""
    global EXCLUSIONS, EXCLUSION_AUDIT
    rows = []
    with gzip.open(ANNOTATION, 'rt') as f:
        for line in f:
            a = line.rstrip().split('\t')
            if a[2].startswith('chr') and a[2][3:].isdigit():
                rows.append(dict(gene=a[12], transcript=a[1], chrom=int(a[2][3:]),
                                 start=int(a[4])+1, end=int(a[5])))
    transcripts = pd.DataFrame(rows)
    EXCLUSIONS = {}; EXCLUSION_AUDIT = {}
    OUT.mkdir(parents=True, exist_ok=True)
    for locus, (chrom, oldlo, oldhi) in LOCI.items():
        t = transcripts[transcripts.chrom.eq(chrom)]
        focal = t[t.gene.isin(FOCAL_GENES[locus])]
        assert set(FOCAL_GENES[locus]) <= set(focal.gene), 'Missing focal annotation'
        centre = (oldlo+oldhi)//2
        corelo = min(oldlo+25000 if locus == 'AMY' else centre, int(focal.start.min()))
        corehi = max(oldhi-25000 if locus == 'AMY' else centre, int(focal.end.max()))
        corelo = (corelo//1000)*1000; corehi = ((corehi+999)//1000)*1000
        while True:
            overlap = t[t.start.le(corehi) & t.end.ge(corelo)]
            # Include EVERY transcript of every gene touching the core.
            genes = t[t.gene.isin(overlap.gene)]
            lo = (min(corelo, int(genes.start.min()))//1000)*1000
            hi = ((max(corehi, int(genes.end.max()))+999)//1000)*1000
            if (lo, hi) == (corelo, corehi): break
            corelo, corehi = lo, hi
        corelo = (corelo//1000)*1000; corehi = ((corehi+999)//1000)*1000
        genes = t[t.gene.isin(t.loc[t.start.le(corehi) & t.end.ge(corelo), 'gene'])]
        assert (genes.start.ge(corelo) & genes.end.le(corehi)).all(), 'Incomplete core after rounding'
        lo, hi = corelo-25000, corehi+25000
        assert (genes.start.ge(lo) & genes.end.le(hi)).all()
        assert lo <= centre <= hi
        EXCLUSIONS[locus] = (chrom, lo, hi)
        EXCLUSION_AUDIT[locus] = dict(original_exclusion=[oldlo, oldhi], corrected_core=[corelo, corehi],
                                      corrected_exclusion=[lo, hi], LD_buffer_bp=25000,
                                      original_focal_transcripts_uncovered=int((~(focal.start.ge(oldlo) & focal.end.le(oldhi))).sum()),
                                      core_genes=sorted(genes.gene.unique()), core_transcripts=len(genes),
                                      uncovered_core_transcripts=0, assertions_passed=True)
        genes.to_csv(OUT/f'{locus}_core_transcripts.tsv', sep='\t', index=False)
        for width in (200000, 1000000):
            outerlo = lo-width//2 if locus == 'AMY' else centre-width//2
            outerhi = hi+width//2 if locus == 'AMY' else centre+width//2
            assert outerlo < centre < outerhi and outerlo < outerhi
            window_genes = t.loc[t.start.le(outerhi) & t.end.ge(outerlo), 'gene']
            inventory = t[t.gene.isin(window_genes)].copy()
            inventory['fully_excluded'] = inventory.start.ge(lo) & inventory.end.le(hi)
            inventory['core_gene'] = inventory.gene.isin(genes.gene)
            assert inventory.loc[inventory.core_gene, 'fully_excluded'].all()
            inventory.to_csv(OUT/f'{locus}_{width}_window_gene_transcripts.tsv', sep='\t', index=False)
            gene_table = inventory.groupby('gene').agg(start=('start', 'min'), end=('end', 'max'),
                           transcripts=('transcript', 'size'), core_gene=('core_gene', 'all'),
                           fully_excluded=('fully_excluded', 'all')).reset_index()
            gene_table.to_csv(OUT/f'{locus}_{width}_window_genes.tsv', sep='\t', index=False)
    save_json(OUT/'annotation_interval_audit.json', EXCLUSION_AUDIT)


def clean_json(value):
    if isinstance(value, dict): return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [clean_json(v) for v in value]
    if isinstance(value, np.ndarray): return clean_json(value.tolist())
    if isinstance(value, np.generic): return clean_json(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def save_json(path, obj):
    path.write_text(json.dumps(clean_json(obj), indent=2, allow_nan=False)+'\n')


def bounds(locus, width):
    if EXCLUSIONS is None: initialize_intervals()
    chrom, lo, hi = EXCLUSIONS[locus]
    centre = (LOCI[locus][1]+LOCI[locus][2])//2
    outerlo = lo-width//2 if locus == 'AMY' else centre-width//2
    outerhi = hi+width//2 if locus == 'AMY' else centre+width//2
    assert outerlo < centre < outerhi and lo <= centre <= hi
    assert EXCLUSION_AUDIT[locus]['uncovered_core_transcripts'] == 0
    return chrom, outerlo, outerhi, lo, hi


def excluded_positions(pos, lo, hi):
    return (np.asarray(pos) >= lo) & (np.asarray(pos) <= hi)


def projection_weights(indices):
    # Per-index seed makes coefficients independent of batch/window membership.
    return np.stack([np.random.default_rng(SEED+int(i)).normal(0, 1/np.sqrt(DIM), DIM)
                     for i in indices]).astype(np.float32)


def encode_block(g, indices):
    """Variant-major ALT dosages; return sketch contribution and marker eligibility."""
    g = np.asarray(g, dtype=np.float32).copy()
    good = np.isfinite(g) & (g >= 0) & (g <= 2)
    counts = good.sum(axis=1)
    means = np.where(good, g, 0).sum(axis=1)/np.maximum(counts, 1)
    residual = np.where(good, g-means[:, None], 0)
    sd = np.sqrt((residual**2).sum(axis=1)/np.maximum(counts, 1))
    eligible = (counts > 0) & (sd > 1e-8)
    if not eligible.any():
        return np.zeros((g.shape[1], DIM), np.float32), eligible, good
    x = residual[eligible]/np.maximum(sd[eligible, None], .1)
    contribution = x.T.astype(np.float32)@projection_weights(np.asarray(indices)[eligible])
    return contribution, eligible, good


def build_sketch(geno, panel, locus):
    """Read all outer-window variants; retain all eligible SNPs in two nested sketches."""
    n = geno.n_samples
    sketches = {w: np.zeros((n, DIM), np.float32) for w in (200000, 1000000)}
    included = np.zeros((n, DIM), np.float32)
    calls = {w: np.zeros(n, np.int32) for w in sketches}
    features = {w: [] for w in sketches}
    total_included = 0
    records = []
    p = panel.reset_index(drop=True)
    for start in range(0, len(p), 256):
        batch = p.iloc[start:start+256]
        indices = batch['index'].to_numpy(np.uint32)
        g = geno.read_variants(indices)
        # Only biallelic SNPs enter the encoding, irrespective of outcomes/labels.
        snp = batch.ref.str.fullmatch('[ACGT]') & batch.alt.str.fullmatch('[ACGT]')
        contribution, eligible, good = encode_block(g, indices)
        eligible &= snp.to_numpy()
        pos = batch.pos.to_numpy()
        _, _, _, lo, hi = bounds(locus, 1000000)
        forbidden = excluded_positions(pos, lo, hi)
        for j, row in enumerate(batch.itertuples(index=False)):
            records.append(dict(index=row.index, chrom=row.chrom, pos=row.pos, ID=row.ID,
                                ref=row.ref, alt=row.alt, excluded=bool(forbidden[j]),
                                eligible_SNP=bool(eligible[j])))
        # Re-encode selected rows so indels and forbidden sites contribute exactly zero.
        if eligible.any():
            c, _, _ = encode_block(g[eligible], indices[eligible])
            included += c
            total_included += int(eligible.sum())
        for w in sketches:
            _, outerlo, outerhi, _, _ = bounds(locus, w)
            keep = eligible & ~forbidden & (pos >= outerlo) & (pos <= outerhi)
            if keep.any():
                c, _, _ = encode_block(g[keep], indices[keep])
                sketches[w] += c
                calls[w] += good[keep].sum(axis=0).astype(np.int32)
                features[w].extend(indices[keep].tolist())
    audit = pd.DataFrame(records)
    audit.to_csv(OUT/f'{locus}_window_variants.tsv', sep='\t', index=False)
    for w in sketches:
        if not features[w]: raise ValueError(f'No usable background features: {locus} {w}')
        sketches[w] /= np.sqrt(len(features[w]))
        calls[w] = calls[w]/len(features[w])
        f = audit[audit['index'].isin(features[w])]
        assert not f.excluded.any()
        f.to_csv(OUT/f'{locus}_{w}_features.tsv', sep='\t', index=False)
        np.save(OUT/f'{locus}_{w}_sketch.npy', sketches[w])
    included /= np.sqrt(max(total_included, 1))
    return sketches, calls, features, included, audit


def fit_clusters(x, k, seed=SEED):
    model = MiniBatchKMeans(n_clusters=k, random_state=seed, n_init=5,
                           batch_size=2048, max_iter=100)
    model.fit(x)
    return model.predict(x)


def choose_k(scores):
    best = max(scores.values())
    return min(k for k, score in scores.items() if score >= best-.01)


def eta_squared(y, labels):
    y = np.asarray(y, float)
    total = ((y-y.mean())**2).sum()
    if total == 0: return 0.
    return float(sum((labels == c).sum()*(y[labels == c].mean()-y.mean())**2
                     for c in np.unique(labels))/total)


def permutation_control(labels):
    rng = np.random.default_rng(SEED+333)
    y = rng.permutation(np.arange(len(labels)) % 2)
    observed = eta_squared(y, labels)
    null = np.array([eta_squared(rng.permutation(y), labels) for _ in range(199)])
    return dict(eta_squared=observed, permutations=199,
                p=(1+int((null >= observed).sum()))/200)


def strip_id(s):
    return re.sub(r'\.(SG|DG|AG|HO|TW|BY|AA|EC|WGC)$', '', str(s))


def panel_metadata(ids):
    a = pd.read_csv(ANNO, sep='\t', dtype=str, keep_default_na=False)
    # Name columns by schema position to preserve published fields without huge headers.
    a = a.iloc[:, [0, 2, 6, 9, 10, 14, 16, 17, 18, 21]].copy()
    a.columns = ['genetic_id', 'biological_id', 'publication', 'date_method', 'BP',
                 'group', 'country', 'latitude', 'longitude', 'data_type']
    direct = {key: part.iloc[0] for key, part in a.groupby('genetic_id') if len(part) == 1}
    aliases = {}
    for _, row in a.iterrows():
        for key in {row.biological_id, strip_id(row.genetic_id)}:
            aliases.setdefault(key, []).append(row)
    rows = []
    for i, sid in enumerate(ids):
        row = direct.get(sid)
        method = 'exact' if row is not None else 'absent'
        if row is None:
            candidates = aliases.get(strip_id(sid), [])
            signatures = {(r.group, r.BP, r.country, r.date_method) for r in candidates}
            if len(signatures) == 1:
                row = candidates[0]; method = 'unique_metadata_signature'
        d = dict(index=i, sample=sid, mapping=method)
        for col in a.columns:
            d[col] = row[col] if row is not None else ''
        d['technology'] = 'shotgun_diploid' if sid.endswith(('.SG', '.DG')) else 'other'
        rows.append(d)
    m = pd.DataFrame(rows)
    for c in ['BP', 'latitude', 'longitude']: m[c] = pd.to_numeric(m[c], errors='coerce')
    m.loc[m.date_method.str.startswith('Modern'), 'BP'] = np.nan
    m.to_csv(OUT/'panel_metadata.tsv', sep='\t', index=False)
    return m


def structural_inventory(ids):
    scratch = OUT/'scratch'; scratch.mkdir(parents=True, exist_ok=True)
    wanted = ['haplotype_structure_for_supplement.tsv', 'ancient_metadata_for_supplement.tsv',
              'haplotype_all_structures.tsv', 'haplotype_struc_dists.tsv',
              'world_wide_cluster_freqs.tsv', 'clusters_processed.tsv',
              'haplotype_representative_structures_rep.tsv', 'haplotype_assignments.txt']
    tables = {}
    with zipfile.ZipFile(ARCHIVE) as z:
        for suffix in wanted:
            hits = [n for n in z.namelist() if n.endswith('/'+suffix)]
            assert len(hits) == 1, (suffix, hits)
            path = scratch/suffix; path.write_bytes(z.read(hits[0]))
            tables[suffix] = pd.read_csv(path, sep='\t')
    h = tables[wanted[0]]
    meta = tables[wanted[1]]
    # Map every archive row to any bundle counterpart, not just selected ancients.
    aliases = {}
    for i, sid in enumerate(ids): aliases.setdefault(strip_id(sid), []).append(i)
    h = h.copy()
    h['panel_index'] = [aliases.get(str(s), [np.nan])[0] if len(aliases.get(str(s), [])) == 1
                        else np.nan for s in h['sample']]
    h = h.merge(meta[['sample', 'age']], on='sample', how='left', validate='many_to_one')
    h['age'] = pd.to_numeric(h.age, errors='coerce')
    h.to_csv(OUT/'structural_mapping.tsv', sep='\t', index=False)
    slots = h.melt(id_vars=['sample', 'type', 'source', 'p1', 'p2', 'panel_index'],
                   value_vars=['H1', 'H2'], var_name='slot', value_name='label').dropna(subset=['label'])
    inventory = slots.groupby(['type', 'source', 'p1', 'p2', 'label'], dropna=False).size().rename('slots').reset_index()
    inventory['frequency_in_partition'] = inventory.slots/inventory.groupby(['type', 'source', 'p1', 'p2'], dropna=False).slots.transform('sum')
    inventory.to_csv(OUT/'structural_partition_frequencies.tsv', sep='\t', index=False)
    counts = slots.groupby('label').size().rename('slots').reset_index()
    counts['frequency'] = counts.slots/len(slots)
    counts['bundle_overlap_slots'] = counts.label.map(slots[slots.panel_index.notna()].groupby('label').size()).fillna(0).astype(int)
    counts.to_csv(OUT/'structural_label_inventory.tsv', sep='\t', index=False)
    structs = tables['haplotype_all_structures.tsv']
    dists = tables['haplotype_struc_dists.tsv']
    # Archive representative-distance strings use a different block nomenclature from
    # all_structures; use the representative table as the EXACT join to distance data.
    reps = tables['haplotype_representative_structures_rep.tsv']
    clusters = tables['clusters_processed.tsv']
    names = tables['haplotype_assignments.txt']
    joined = structs.merge(clusters, left_on='contig', right_on='sample', how='left')
    joined = joined.merge(names[['chrom', 'sample_id', 'haplotype']], left_on='contig', right_on='chrom',
                          how='left', validate='one_to_one')
    # Representative contig intervals have a few boundary updates; the stable
    # contig identifier before ':' is unique in BOTH tables and carries no label inference.
    names = names.copy(); reps = reps.copy()
    names['stable_contig'] = names.chrom.str.split(':').str[0]
    reps['stable_contig'] = reps.contig.str.split(':').str[0]
    assert names.stable_contig.is_unique and reps.stable_contig.is_unique
    rep_labels = reps.merge(names[['stable_contig', 'haplotype']], on='stable_contig', how='left', validate='one_to_one')
    rep_labels.to_csv(OUT/'distance_representative_labels.tsv', sep='\t', index=False)
    join_identity = joined.hap_struc_d.eq(names.set_index('chrom').loc[joined.contig, 'hap_struc_d'].to_numpy())
    assert join_identity.all(), 'Contig-to-label structure strings must match exactly'
    region_lookup = h.drop_duplicates('sample').set_index('sample').p1.to_dict()
    joined['population_partition'] = [region_lookup.get(re.sub(r'^AV_', '', str(s)), 'UNKNOWN') for s in joined.sample_id]
    joined.to_csv(OUT/'sequence_structure_inventory.tsv', sep='\t', index=False)
    label_rows = []
    for label, part in joined.groupby('haplotype'):
        rr = rep_labels[rep_labels.haplotype.eq(label)]
        within = dists[dists.hap1.isin(rr.hap_struc_d) & dists.hap2.isin(rr.hap_struc_d)]
        label_rows.append(dict(label=label, assembled_contigs=len(part),
                               exact_structure_strings=part.hap_struc_d.nunique(),
                               known_population_partitions=part.loc[part.population_partition.ne('UNKNOWN'), 'population_partition'].nunique(),
                               distance_representatives=len(rr), within_label_distance_pairs=len(within),
                               within_label_min_edit_d=float(within.edit_d.min()) if len(within) else None,
                               within_label_max_edit_d=float(within.edit_d.max()) if len(within) else None))
    pd.DataFrame(label_rows).to_csv(OUT/'structural_label_identity_inventory.tsv', sep='\t', index=False)
    seq_summary = dict(contigs=len(structs), unique_structures=structs.hap_struc_d.nunique(),
                       distance_rows=len(dists), zero_distance_rows=int((dists.edit_d == 0).sum()),
                       positive_distance_rows=int((dists.edit_d > 0).sum()),
                       max_edit_distance=float(dists.edit_d.max()),
                       numeric_cluster_mapped_contigs=int(joined.cluster.notna().sum()),
                       representatives=len(reps),
                       representative_strings_in_distance=int(reps.hap_struc_d.isin(set(dists.hap1)|set(dists.hap2)).sum()),
                       label_mapped_contigs=int(joined.haplotype.notna().sum()),
                       label_mapped_distance_representatives=int(rep_labels.haplotype.notna().sum()),
                       label_distance_join='haplotype_assignments.txt exact contig/string join; representative stable-contig IDs',
                       origin_recurrence='not identified: structure heterogeneity is not independent molecular origin')
    return h, dict(archive_individual_rows=len(h), labelled_slots=len(slots), labels=len(counts),
                   bundle_matched_individual_rows=int(h.panel_index.notna().sum()),
                   bundle_labelled_slots=int(slots.panel_index.notna().sum()),
                   dated_ancient_bundle_rows=int((h.type.eq('Ancient') & h.age.notna() & h.panel_index.notna()).sum()),
                   structure_distances=seq_summary)


def composition(m, labels, locus):
    x = m.copy(); x['cluster'] = labels
    for field in ['group', 'country', 'publication', 'technology', 'data_type']:
        x[field] = x[field].replace('', 'UNMAPPED')
        x.groupby(['cluster', field]).size().rename('n').reset_index().to_csv(
            OUT/f'{locus}_composition_{field}.tsv', sep='\t', index=False)
    summary = x.groupby('cluster').agg(n=('sample', 'size'), mapped=('mapping', lambda s: (s != 'absent').sum()),
                                       dated_n=('BP', 'count'), min_BP=('BP', 'min'), median_BP=('BP', 'median'),
                                       max_BP=('BP', 'max'), median_latitude=('latitude', 'median'),
                                       median_longitude=('longitude', 'median'))
    summary.to_csv(OUT/f'{locus}_composition_summary.tsv', sep='\t')
    known = m.group.ne('') & m.group.ne('..')
    country = m.country.ne('') & m.country.ne('..')
    batch_values, batch_ns = [], []
    for _, part in x[known].groupby('group'):
        if len(part) >= 20 and part.technology.nunique() >= 2:
            batch_values.append(normalized_mutual_info_score(part.cluster, part.technology))
            batch_ns.append(len(part))
    return dict(mapped_group_n=int(known.sum()), mapped_country_n=int(country.sum()),
                group_NMI=normalized_mutual_info_score(m.loc[known, 'group'], labels[known]),
                country_NMI=normalized_mutual_info_score(m.loc[country, 'country'], labels[country]),
                within_group_technology_NMI=float(np.average(batch_values, weights=batch_ns)) if batch_ns else None,
                within_group_technology_n=sum(batch_ns), within_group_technology_groups=len(batch_ns))


def design_floor(time_axis, labels, p):
    """Outcome-blind independent-binomial floor conditional on cluster intercepts."""
    t = np.asarray(time_axis, float); labels = np.asarray(labels)
    if len(t) < 30 or np.ptp(t) < 2 or not 0 < p < 1:
        return dict(n=len(t), frequency_assumption=p, identifiable=False)
    tr = t.copy()
    for c in np.unique(labels): tr[labels == c] -= t[labels == c].mean()
    information = float(tr@tr)
    if information <= 0: return dict(n=len(t), identifiable=False)
    floor = 2.802*np.sqrt(p*(1-p)/(2*information))
    return dict(n=len(t), frequency_assumption=p, time_information_kyr2=information,
                age_span_kyr=float(np.ptp(t)), identifiable=True,
                MDE80_ideal_per_kyr=floor, MDE80_ideal_2kyr=2*floor)


def temporal_fit(y, time_axis, labels):
    y = np.asarray(y, float); t = np.asarray(time_axis, float); labels = np.asarray(labels)
    ideal = design_floor(t, labels, float(y.mean()))
    if not ideal['identifiable']: return dict(**ideal, HC3_identifiable=False)
    categories = np.unique(labels)
    x = np.column_stack([(labels == c).astype(float) for c in categories]+[t])
    if np.linalg.matrix_rank(x) < x.shape[1] or len(y) <= x.shape[1]+2:
        return dict(**ideal, HC3_identifiable=False)
    inv = np.linalg.inv(x.T@x); beta = inv@x.T@y
    h = np.sum((x@inv)*x, axis=1)
    if (h >= 1-1e-8).any(): return dict(**ideal, HC3_identifiable=False)
    r = (y-x@beta)/(1-h)
    cov = inv@(x.T@(r[:, None]**2*x))@inv
    se = np.sqrt(max(cov[-1, -1], 0)); df = len(y)-x.shape[1]
    floor = (student_t.ppf(.975, df)+norm.ppf(.8))*se
    return dict(**ideal, HC3_identifiable=True, MDE80_HC3_per_kyr=floor, MDE80_HC3_2kyr=2*floor,
                beta_per_kyr=beta[-1], HC3_SE=se,
                p=2*student_t.sf(abs(beta[-1]/se), df) if se > 0 else None)


def outcome_tests(geno, panel, locus, labels, m, structural):
    outcomes = {}
    if locus == 'AMY':
        h = structural[structural.panel_index.notna()].copy()
        for label in sorted(set(h.H1.dropna())|set(h.H2.dropna())):
            y = (h.H1.eq(label).astype(int)+h.H2.eq(label).astype(int))/2
            outcomes[label] = (h.panel_index.astype(int).to_numpy(), y.to_numpy(), h.age.to_numpy(float))
        copies1 = h.H1.str.extract(r'^H(\d+)')[0].astype(float)
        copies2 = h.H2.str.extract(r'^H(\d+)')[0].astype(float)
        for cut in [1, 2]:
            y = ((copies1 > cut).astype(int)+(copies2 > cut).astype(int))/2
            outcomes[f'B_A_le{cut}'] = (h.panel_index.astype(int).to_numpy(), y.to_numpy(), h.age.to_numpy(float))
    else:
        centre = (LOCI[locus][1]+LOCI[locus][2])//2
        hit = panel[(panel.pos == centre) & panel.ref.str.fullmatch('[ACGT]') & panel.alt.str.fullmatch('[ACGT]')]
        if len(hit) != 1: return dict(status='tested SNP unavailable/ambiguous', matches=len(hit))
        g = geno.read_variants(hit['index'].to_numpy(np.uint32))[0]
        usable = np.isfinite(g) & (g >= 0) & (g <= 2)
        ix = np.flatnonzero(usable)
        outcomes[str(hit.iloc[0].ID)] = (ix, g[ix]/2, m.BP.to_numpy(float)[ix])
    frequencies, temporal, sharing = [], [], []
    for name, (ix, y, age) in outcomes.items():
        labs = labels[ix]; t = -age/1000; dated = np.isfinite(t)
        row = temporal_fit(y[dated], t[dated], labs[dated]) if dated.any() else dict(HC3_identifiable=False, n=0)
        temporal.append(dict(variant=name, scope='common_within', **row))
        allele_counts = []
        supported = 0
        for c in np.unique(labs):
            mask = labs == c; nn = int(mask.sum()); copies = float(2*y[mask].sum())
            allele_counts.append(copies)
            supported += int(copies >= 5 and 2*nn-copies >= 5)
            frequencies.append(dict(variant=name, cluster=int(c), n=nn, slots=2*nn,
                                    expected_ALT_copies=copies, frequency=float(y[mask].mean()),
                                    dated_n=int((mask & dated).sum())))
            tt = mask & dated
            row = temporal_fit(y[tt], t[tt], np.zeros(tt.sum(), int)) if tt.any() else dict(HC3_identifiable=False, n=0)
            temporal.append(dict(variant=name, scope=f'cluster_{c}', **row))
        ac = np.asarray(allele_counts); shares = ac/ac.sum() if ac.sum() else ac
        effective = float(np.exp(-np.sum(shares[shares > 0]*np.log(shares[shares > 0])))) if ac.sum() else 0
        sharing.append(dict(variant=name, supported_backgrounds=supported, effective_backgrounds=effective,
                            independent_origin_inference='not identifiable from unphased cluster sharing'))
    pd.DataFrame(frequencies).to_csv(OUT/f'{locus}_frequencies.tsv', sep='\t', index=False)
    pd.DataFrame(temporal).to_csv(OUT/f'{locus}_temporal.tsv', sep='\t', index=False)
    pd.DataFrame(sharing).to_csv(OUT/f'{locus}_background_sharing.tsv', sep='\t', index=False)
    return dict(status='diagnostic outcomes completed', outcomes=len(outcomes), temporal=temporal, sharing=sharing)


def analyze_locus(geno, panel, locus, m, structural):
    started = time.monotonic()
    sketches, calls, features, included, audit = build_sketch(geno, panel, locus)
    rng = np.random.default_rng(SEED)
    evaluation = np.sort(rng.choice(geno.n_samples, min(1500, geno.n_samples), replace=False))
    assignments = {}; scores = {}; chosen = {}
    for w, x in sketches.items():
        scores[w] = {}
        for k in KS:
            labels = fit_clusters(x, k)
            assignments[w, k] = labels
            score = silhouette_score(x[evaluation], labels[evaluation]) if len(np.unique(labels[evaluation])) > 1 else -1
            scores[w][k] = float(score)
        chosen[w] = choose_k(scores[w])
    k = chosen[1000000]; labels = assignments[1000000, k]
    neighbor = [q for q in KS if abs(q-k) == 2]
    stability = dict(same_K_window_ARI=adjusted_rand_score(labels, assignments[200000, k]),
                     selected_K_window_ARI=adjusted_rand_score(labels, assignments[200000, chosen[200000]]),
                     neighbor_K_ARI={q: adjusted_rand_score(labels, assignments[1000000, q]) for q in neighbor},
                     second_seed_ARI=adjusted_rand_score(labels, fit_clusters(sketches[1000000], k, SEED+1)),
                     exclusion_off_ARI=adjusted_rand_score(labels, fit_clusters(included, k)))
    comp = composition(m, labels, locus)
    control = permutation_control(labels)
    coverage = {w: dict(features=len(features[w]), assigned=geno.n_samples,
                        observed_95pct_n=int((calls[w] >= .95).sum()),
                        observed_95pct_fraction=float((calls[w] >= .95).mean()),
                        min_call_fraction=float(calls[w].min()), median_call_fraction=float(np.median(calls[w]))) for w in sketches}
    forbidden_in_features = int(audit[audit['index'].isin(features[1000000])].excluded.sum())
    exclusion = dict(excluded_window_variants=int(audit.excluded.sum()), forbidden_primary_features=forbidden_in_features,
                     primary_feature_count=len(features[1000000]),
                     feature_index_sha256=hashlib.sha256(np.asarray(features[1000000], np.uint32).tobytes()).hexdigest())
    gates = dict(a=bool(scores[1000000][k] >= .1 and comp['group_NMI'] >= .1 and comp['country_NMI'] >= .05
                        and comp['within_group_technology_NMI'] is not None and comp['within_group_technology_NMI'] <= .2),
                 b=bool(stability['same_K_window_ARI'] >= .65 and stability['second_seed_ARI'] >= .8
                        and all(v >= .6 for v in stability['neighbor_K_ARI'].values())),
                 c=bool(forbidden_in_features == 0 and control['p'] >= .01),
                 d=bool(all(v['observed_95pct_fraction'] >= .99 for v in coverage.values())))
    assignment_frame = m[['index', 'sample']].copy()
    for (w, kk), ll in assignments.items(): assignment_frame[f'window_{w}_K{kk}'] = ll
    for w, cf in calls.items(): assignment_frame[f'call_fraction_{w}'] = cf
    assignment_frame.to_csv(OUT/f'{locus}_assignments.tsv', sep='\t', index=False)
    result = dict(locus=locus, interval_audit=EXCLUSION_AUDIT[locus], scores=scores, chosen_K=chosen, composition=comp,
                  stability=stability, negative_control=control, exclusion=exclusion,
                  coverage=coverage, gates=gates, all_gates_pass=all(gates.values()))
    if locus == 'AMY':
        hh = structural[structural.type.eq('Ancient') & structural.panel_index.notna() & structural.age.notna()]
        idx = hh.panel_index.astype(int).to_numpy()
        result['outcome_blind_design_floor'] = {
            'B_A_le1_published_p_0.842222': design_floor(-hh.age.to_numpy(float)/1000, labels[idx], .842222),
            'B_A_le2_published_p_0.742222': design_floor(-hh.age.to_numpy(float)/1000, labels[idx], .742222)}
    if all(gates.values()):
        result['outcomes'] = outcome_tests(geno, panel, locus, labels, m, structural)
    else:
        result['outcomes'] = dict(status='STOPPED: validation failed', failed_gates=[g for g, ok in gates.items() if not ok])
    result['elapsed_seconds'] = time.monotonic()-started
    result['peak_RSS_MiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    save_json(OUT/f'{locus}_summary.json', result)
    print(locus, json.dumps(clean_json(dict(chosen_K=k, gates=gates, stability=stability,
                                          elapsed_seconds=result['elapsed_seconds']))), flush=True)
    return result


def scan_windows():
    records = {locus: [] for locus in LOCI}
    # Stream full immutable pvar once. Cache only this task's focal window coordinates.
    windows = {bounds(locus, 1000000)[0]: (locus, bounds(locus, 1000000)[1],
                                         bounds(locus, 1000000)[2]) for locus in LOCI}
    with open(paths()[1]) as f:
        index = 0
        for line in f:
            if line.startswith('#'): continue
            chrom, pos, vid, ref, alt = line.rstrip('\n').split('\t')[:5]
            chrom = int(chrom)
            window = windows.get(chrom)
            if window is not None:
                locus, lo, hi = window
                pos = int(pos)
                if lo <= pos <= hi:
                    records[locus].append((index, chrom, pos, vid, ref, alt))
            index += 1
    for locus, values in records.items():
        pd.DataFrame(values, columns=['index', 'chrom', 'pos', 'ID', 'ref', 'alt']).to_csv(
            OUT/f'{locus}_pvar.tsv', sep='\t', index=False)
    return index


def record_provenance(summary):
    from scripts import provenance as p
    p.get_git_info = lambda: {'commit': 'not inspected: no-git operator rule', 'branch': 'not inspected', 'dirty': 'unknown'}
    r = p.RunRecorder('haplotype-background', command=[sys.executable]+sys.argv)
    r.set_seed(SEED); r.record_parameter('summary', summary)
    for f in [PRE, Path(__file__), ANNO, ANNOTATION, ARCHIVE, paths()[2]]: r.record_input(f)
    for f in paths()[:2]:
        st = f.stat(); r.record_parameter(str(f), dict(size_bytes=st.st_size, mtime_ns=st.st_mtime_ns,
                                                       hash='not computed: bounded immutable-input IO'))
    for f in OUT.glob('*'):
        if f.is_file(): r.record_output(f)
    r.save()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--scan-only', action='store_true')
    ap.add_argument('--loci', nargs='+', default=list(LOCI), choices=list(LOCI))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    initialize_intervals()
    start = datetime.now(timezone.utc).isoformat(); timer = time.monotonic()
    summary = dict(started_utc=start, prereg_mtime_ns=PRE.stat().st_mtime_ns,
                   prereg_sha256=hashlib.sha256(PRE.read_bytes()).hexdigest(), status='running', loci={})
    save_json(OUT/'summary.json', summary)
    if not all((OUT/f'{loc}_pvar.tsv').exists() for loc in LOCI):
        summary['pvar_variants_counted'] = scan_windows()
        save_json(OUT/'scan_summary.json', summary)
    if args.scan_only:
        record_provenance(summary); return
    ids = parse_psam(paths()[2])
    structural, inventory = structural_inventory(ids)
    metadata = panel_metadata(ids)
    summary['structural_inventory'] = inventory
    with AkbariGeno() as geno:
        assert geno.n_samples == len(ids)
        summary['panel'] = dict(samples=geno.n_samples, variants=geno.n_variants,
                                nominal_chromosome_slots=2*geno.n_samples,
                                sample_subset='NONE: every psam row', phased_haplotypes='not inferred')
        save_json(OUT/'summary.json', summary)
        for locus in args.loci:
            panel = pd.read_csv(OUT/f'{locus}_pvar.tsv', sep='\t')
            summary['loci'][locus] = analyze_locus(geno, panel, locus, metadata, structural)
            save_json(OUT/'summary.json', summary)
    summary['status'] = 'complete; failed loci stopped individually'
    summary['elapsed_seconds'] = time.monotonic()-timer
    save_json(OUT/'summary.json', summary)
    record_provenance(summary)


if __name__ == '__main__': main()
