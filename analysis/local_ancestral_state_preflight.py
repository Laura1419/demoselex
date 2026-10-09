#!/usr/bin/env python3
"""Fail-closed input/evidence audit; does not infer states or fit AMY outcomes."""
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
from analysis.local_ancestral_state_checks import ideal_frequency_floor, spectral_effective_count
from scripts import provenance

# Operator prohibition: disable git BEFORE RunRecorder constructor.
provenance.get_git_info = lambda: {"disabled": True, "reason": "operator forbids all git commands"}
OUT = Path('data/derived/local_ancestral_state')
PRE = Path('docs/local-ancestral-state-prereg.md')
METHOD = Path('docs/local-ancestral-state-method.md')
PREFIX = Path('data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes')


def main():
    assert PRE.exists() and METHOD.exists(), 'Study and registration must exist first'
    OUT.mkdir(parents=True, exist_ok=True)
    prereg = PRE.read_bytes()
    digest = hashlib.sha256(prereg).hexdigest()
    snapshot = OUT / 'prereg.before_preflight.md'
    if snapshot.exists():
        assert snapshot.read_bytes() == prereg, 'Prereg changed: append and explicitly register amendment'
    else:
        snapshot.write_bytes(prereg)
    registration = dict(prereg_sha256=digest, snapshot=str(snapshot),
                        snapshotted_utc=datetime.now(timezone.utc).isoformat(),
                        real_data_outcome_statistics_run=False)
    ledger = OUT / 'registration.json'
    if not ledger.exists():
        ledger.write_text(json.dumps(registration, indent=2)+'\n')
    else:
        assert json.loads(ledger.read_text())['prereg_sha256'] == digest
    recorder = provenance.RunRecorder('local-ancestral-state-preflight')
    for path in [__file__, PRE, METHOD, snapshot, ledger,
                 'analysis/local_ancestral_state_checks.py',
                 'analysis/test_local_ancestral_state_checks.py']:
        recorder.record_input(path)

    def load_json(path):
        recorder.record_input(path)
        return json.loads(Path(path).read_text())

    # Metadata counts only. Do not read any genotype/variant outcome.
    psam = Path(str(PREFIX)+'.psam'); pgen = Path(str(PREFIX)+'.pgen')
    recorder.record_input(psam)
    sample_ids = [line.split('\t')[0] for line in psam.read_text().splitlines()
                  if line and not line.startswith('#')]
    reader = pgenlib.PgenReader(str(pgen).encode())
    n, variants = reader.get_raw_sample_ct(), reader.get_variant_ct()
    reader.close()
    assert n == len(sample_ids) == 16339 and variants == 52382872
    recorder.record_parameter('large_input_size_mtime_only', {
        str(pgen): {'size': pgen.stat().st_size, 'mtime_ns': pgen.stat().st_mtime_ns}})

    metadata_path = Path('data/derived/haplotype_background/panel_metadata.tsv')
    recorder.record_input(metadata_path)
    meta = pd.read_csv(metadata_path, sep='\t')
    assert len(meta) == n and meta['sample'].tolist() == sample_ids
    mapped = meta.genetic_id.notna()
    date_valid = pd.to_numeric(meta.BP, errors='coerce').gt(0)

    annotation = Path('data/raw/annotation/hg19_refGene.txt.gz')
    recorder.record_input(annotation)
    amy_rows = []
    with gzip.open(annotation, 'rt') as stream:
        for line in stream:
            fields = line.rstrip().split('\t')
            if fields[2] == 'chr1' and fields[12] in ['AMY1A','AMY1B','AMY1C','AMY2A','AMY2B']:
                amy_rows.append((fields[12], int(fields[4])+1, int(fields[5])))
    assert amy_rows and max(row[2] for row in amy_rows if row[0].startswith('AMY1')) == 104301314
    uncovered = sum(not (103915000 <= row[1] <= row[2] <= 104327000) for row in amy_rows)
    assert uncovered == 0

    # Saved precedent evidence is explicitly not a newly fitted validation.
    prior = load_json('data/derived/haplotype_background/AMY_summary.json')
    local = load_json('data/derived/local_ancestry/validation.json')
    confusion_path = Path('data/derived/local_ancestry/confusion.tsv')
    recorder.record_input(confusion_path)
    confusion = pd.read_csv(confusion_path, sep='\t').fillna(0)
    global_baseline = load_json('data/derived/local_ancestry/verified_global_baseline.json')
    released = load_json('data/derived/local_ancestry/verified_akbari_windows.json')
    means_path = Path('data/derived/local_ancestry/verified_CN_means.tsv')
    recorder.record_input(means_path)
    cn_means = pd.read_csv(means_path, sep='\t').to_dict('records')
    relatedness_path = Path('data/derived/two_channel/amy_relatedness.tsv')
    recorder.record_input(relatedness_path)
    old_r = pd.read_csv(relatedness_path, sep='\t', index_col=0)
    # The old seven-component calculation excludes the single CHG carrier.
    labels = [label for label in old_r.index if label != 'CHG']
    r = old_r.loc[labels, labels].to_numpy()
    old_summary = load_json('data/derived/two_channel/amy_summary.json')
    old_neff = float(len(labels)**2 / r.sum())
    assert np.isclose(old_neff, old_summary['ordering']['ordering_neff'])

    results = {
        'status': 'STOPPED_BEFORE_STATE_INFERENCE',
        'registration': json.loads(ledger.read_text()),
        'new_biological_statistic': None,
        'new_assignment_confusion_matrix': None,
        'new_assignment_stability': None,
        'new_assignment_usable_n': 0,
        'validated_evaluable_shared_states': 0,
        'biological_shared_state_count': None,
        'local_shared_state_N_eff': None,
        'achieved_local_detection_floor': None,
        'input_metadata': {'panel_n': n, 'variants': variants, 'AADR_mapped_n': int(mapped.sum()),
                           'AADR_date_positive_n': int(date_valid.sum()),
                           'AADR_country_mapped_n': int(meta.country.notna().sum()),
                           'AMY1_annotation_end': 104301314, 'uncovered_focal_transcripts': uncovered},
        'standard_tool_in_PATH': {name: shutil.which(name) for name in ['elai','elai-lin','rfmix','flare']},
        'blockers': ['No validated lower-state joint posterior/central-gap export',
                     'No calibrated genotype uncertainty/QC or masked-preprocessing leakage evidence',
                     'No local truth, neutral calibration or shared-state power validation'],
        'precedent_NOT_new_validation': {
            'attempt1': {'validation': local, 'confusion': confusion.to_dict('records')},
            'attempt2_AMY': prior,
            'saved_global_baseline_NOT_local': global_baseline,
            'saved_released_windows_NOT_new_scan': released,
            'saved_CN_means_NOT_new_local_outcome': cn_means,
            'seven_global_component_N_eff_old_definition': old_neff,
            'same_global_matrix_spectral_count_different_definition': spectral_effective_count(r)},
        'synthetic_design_floors_NOT_achieved_power': {
            'n225_balanced_p_half_Delta': ideal_frequency_floor(112,113),
            'n5737_balanced_p_half_Delta': ideal_frequency_floor(2868,2869),
            'n22_balanced_p_half_Delta': ideal_frequency_floor(11,11)},
        'no_git_commands': True, 'no_files_staged_by_worker': True,
        'repository_index_state': 'not inspected; operator prohibits git commands',
    }
    output = OUT / 'preflight.json'
    output.write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    recorder.record_output(output)
    recorder.record_parameter('read_scope', 'metadata and saved precedent evidence; no genotype outcomes read')
    recorder.save()
    print(json.dumps({key: value for key,value in results.items()
                      if key != 'precedent_NOT_new_validation'}, indent=2))


if __name__ == '__main__':
    main()
