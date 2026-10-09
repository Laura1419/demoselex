#!/usr/bin/env python3
"""Small, streaming continuation audit; never infer states or fit biological outcomes."""
from __future__ import annotations

from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd
import pgenlib

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis.local_ancestral_state_checks import AMY_EXCLUSION, assert_amy_features_excluded
from analysis.load_akbari import iter_chunks
from scripts import provenance

# Disable the repository library's git subprocesses BEFORE constructing any recorder.
provenance.get_git_info = lambda: {"commit": "not-probed", "branch": "not-probed", "dirty": None,
                                    "disabled": True, "reason": "operator prohibits all git commands"}
OUT = Path('data/derived/local_ancestral_state')
PRE = Path('docs/local-ancestral-state-prereg.md')
PREFIX = Path('data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes')
TOOL_REVISION = '4abc4a8fcb824d20db2dde672cfaf5c92fd09853'


def main():
    start = datetime.now(timezone.utc).isoformat()
    free_start = shutil.disk_usage('/').free
    print(f'[observation] Start UTC={start}; root free bytes={free_start}', flush=True)
    original = (OUT / 'prereg.before_preflight.md').read_bytes()
    prereg = PRE.read_bytes()
    assert prereg.startswith(original), 'Original registration must be preserved verbatim'
    assert b'CONTINUATION REGISTRATION' in prereg
    snapshot = OUT / 'prereg.before_resume_audit.md'
    if snapshot.exists():
        assert snapshot.read_bytes() == prereg, 'Append-only amendment needs a new audit registration'
    else:
        snapshot.write_bytes(prereg)
    rec = provenance.RunRecorder('local-ancestral-state-resume-audit')
    rec.set_seed(0)  # Deterministic audit; no random state inference.

    def read_json(path):
        rec.record_input(path)
        return json.loads(Path(path).read_text())

    original_ledger = read_json(OUT / 'registration.json')
    rec.record_input(OUT / 'prereg.before_preflight.md')
    assert hashlib.sha256(original).hexdigest() == original_ledger['prereg_sha256']
    for path in [__file__, PRE, snapshot, 'docs/local-ancestral-state-method.md',
                 'analysis/local_ancestral_state_checks.py',
                 'analysis/test_local_ancestral_state_checks.py', 'analysis/load_akbari.py']:
        rec.record_input(path)
    lo, hi = AMY_EXCLUSION
    assert (lo, hi) == (103915000, 104327000)
    prior = read_json('data/derived/haplotype_background/AMY_summary.json')
    assert prior['interval_audit']['corrected_exclusion'] == [lo, hi]
    core_genes = set(prior['interval_audit']['core_genes'])
    annotation = Path('data/raw/annotation/hg19_refGene.txt.gz')
    rec.record_input(annotation)
    transcripts = []
    with gzip.open(annotation, 'rt') as stream:
        for line in stream:
            f = line.rstrip().split('\t')
            if f[2] == 'chr1' and f[12] in core_genes:
                transcripts.append((f[12], int(f[4])+1, int(f[5])))
    assert len(transcripts) == prior['interval_audit']['core_transcripts'] == 22
    assert all(lo <= start <= end <= hi for _, start, end in transcripts)
    amy1_end = max(end for gene, _, end in transcripts if gene.startswith('AMY1'))
    assert amy1_end == 104301314
    forbidden_path = Path('data/derived/haplotype_background/AMY_window_variants.tsv')
    rec.record_input(forbidden_path)
    variants = pd.read_csv(forbidden_path, sep='\t')
    forbidden = set(variants.loc[variants.chrom.eq(1) & variants.pos.between(lo, hi), 'index'])
    feature_audit = {}
    for width in [200000, 1000000]:
        path = Path(f'data/derived/haplotype_background/AMY_{width}_features.tsv')
        rec.record_input(path)
        features = pd.read_csv(path, sep='\t')
        assert_amy_features_excluded(features.chrom.to_numpy(), features.pos.to_numpy())
        intersections = len(set(features['index']) & forbidden)
        assert intersections == 0
        feature_audit[str(width)] = {'features': len(features), 'forbidden_index_intersections': intersections}

    psam = Path(str(PREFIX)+'.psam')
    pgen = Path(str(PREFIX)+'.pgen')
    rec.record_input(psam)
    reader = pgenlib.PgenReader(str(pgen).encode())
    try:
        n, m = reader.get_raw_sample_ct(), reader.get_variant_ct()
    finally:
        reader.close()
    assert (n, m) == (16339, 52382872)
    sample_ids = [line.split()[0] for line in psam.read_text().splitlines() if not line.startswith('#')]
    meta_path = Path('data/derived/haplotype_background/panel_metadata.tsv')
    rec.record_input(meta_path)
    meta = pd.read_csv(meta_path, sep='\t')
    assert meta['sample'].tolist() == sample_ids and len(meta) == n
    dated = meta[meta.BP.gt(0)]
    cohort_path = Path('data/derived/amy_binarized/primary_cohort.tsv')
    rec.record_input(cohort_path)
    cohort = pd.read_csv(cohort_path, sep='\t')
    assert len(cohort) == cohort['sample'].nunique() == 225
    bundle_set = set(sample_ids)
    assert cohort.pgen_id.isin(bundle_set).all()
    old = read_json('data/derived/amy_binarized/summary.json')
    saved = read_json('data/derived/local_ancestry/verified_global_baseline.json')
    results = {(row['locus'], row['model']): row for row in old['results']}
    primary_on, primary_off = [results[('primary_le1', model)] for model in ['GRM_on', 'GRM_off']]
    retention = abs(primary_on['z_score']/primary_off['z_score'])
    le2_x = abs(results[('sensitivity_le2', 'GRM_on')]['X_equivalent'])
    assert np.isclose(retention, saved['saved_GLMM']['primary_GRM_retention'])
    assert np.isclose(le2_x, saved['saved_GLMM']['sensitivity_le2_X_equivalent'])
    for cut in ['primary_le1', 'sensitivity_le2']:
        for field, saved_field in [('beta_within', 'beta_time'), ('MDE80_HC3_per_kyr', 'MDE80_HC3'),
                                    ('carriage_absolute_share', 'carriage_share')]:
            assert np.isclose(old['channels'][cut][field], saved[cut][saved_field])

    windows = {'AMY': (1, 103500000, 104500000), 'LCT': (2, 136300000, 136900000),
               'FADS': (11, 61300000, 61700000), 'SLC45A2': (5, 33800000, 34200000),
               'SLC24A5': (15, 48200000, 48600000)}
    audit = {name: {'n': 0, 'max_abs_X': 0., 'above_5_45': 0} for name in windows}
    for chunk in iter_chunks(chunksize=100000, usecols=['CHROM', 'POS', 'X']):
        for name, (chrom, lower, upper) in windows.items():
            x = chunk.loc[chunk.CHROM.eq(chrom) & chunk.POS.between(lower, upper), 'X'].abs()
            if len(x):
                audit[name]['n'] += len(x)
                audit[name]['max_abs_X'] = max(audit[name]['max_abs_X'], float(x.max()))
                audit[name]['above_5_45'] += int(x.gt(5.45).sum())
    expected = read_json('data/derived/local_ancestry/verified_akbari_windows.json')
    assert audit == expected, 'Fresh streamed historical-statistic audit disagrees with saved evidence'
    for large in [pgen, Path('data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz')]:
        rec.record_parameter(str(large), {'size': large.stat().st_size, 'mtime_ns': large.stat().st_mtime_ns,
                                       'full_hash_omitted': 'avoid giant file reads; raw input unchanged'})
    source_evidence = []
    for path in sorted((OUT / 'literature').glob('ELAI_4abc4a8f_*.json')):
        source = read_json(path)
        assert source['source_revision'] == TOOL_REVISION
        assert hashlib.sha256(source['body_utf8'].encode()).hexdigest() == source['sha256']
        if 'control.cpp' in path.name:
            assert 'ELAI version 1.01' in source['body_utf8']
        source_evidence.append({'path': str(path), 'url': source['url'], 'sha256': source['sha256']})
    assert len(source_evidence) == 3
    output = OUT / 'resume_audit.json'
    report = {
        'status': 'STOPPED_BEFORE_STATE_INFERENCE', 'started_utc': start,
        'prereg_sha256': hashlib.sha256(prereg).hexdigest(), 'prereg_snapshot': str(snapshot),
        'tool_source_revision': TOOL_REVISION, 'tool_source_version_banner': 'ELAI version 1.01',
        'source_evidence': source_evidence,
        'tool_commit_metadata_url': f'https://api.github.com/repos/haplotype/ELAI/commits/{TOOL_REVISION}',
        'exclusion': {'inclusive_GRCh37': [lo, hi], 'AMY1_end': amy1_end,
                      'core_transcripts': len(transcripts), 'uncovered_core_transcripts': 0,
                      'saved_feature_manifest_audit_NOT_new_assignment': feature_audit},
        'input_counts': {'panel_n': n, 'variants': m, 'dated_metadata_rows': len(dated),
                         'dated_unique_biological_ids': dated.biological_id.nunique(),
                         'dated_3500_6500_BP_rows': int(dated.BP.between(3500,6500).sum()),
                         'measured_structural_dated_bundle_n': len(cohort)},
        'carried_facts': {'fresh_streamed_released_windows': audit,
                          'saved_primary_GRM_z_retention': retention, 'saved_A_le2_abs_X_equivalent': le2_x,
                          'saved_frequency_channels_NOT_new_local_statistics': saved},
        'previous_failed_method_NOT_new_validation': {
            'structure': prior['composition'], 'silhouettes': prior['scores'],
            'stability': prior['stability'], 'coverage': prior['coverage'],
            'negative_control': prior['negative_control'], 'gates': prior['gates']},
        'new_method_gates': {name: 'NOT_RUN_NOT_PASSED' for name in
                             ['implementation', 'structure', 'local_truth', 'stability', 'coverage',
                              'negative_control_leakage', 'identifiability_null', 'sharing_power']},
        'new_biological_statistic': None, 'local_shared_state_N_eff': None,
        'achieved_local_detection_floor': None, 'validated_evaluable_shared_states': 0,
        'root_free_bytes_start': free_start, 'root_free_bytes_end': shutil.disk_usage('/').free,
        'no_git_commands': True, 'no_files_staged_by_worker': True,
        'repository_index_state': 'not inspected: all git commands prohibited',
        'scratch_created': [],
    }
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    rec.record_output(output)
    rec.record_output(snapshot)
    rec.record_parameter('read_scope', 'metadata, saved summaries and chunked historical released X; no genotype arrays')
    rec.save()
    print(f'[observation] End root free bytes={report["root_free_bytes_end"]}; STOPPED; no local statistic', flush=True)
    print(json.dumps({'windows': audit, 'counts': report['input_counts']}, indent=2))


if __name__ == '__main__':
    main()
