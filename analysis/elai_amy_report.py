#!/usr/bin/env python3
"""Report ELAI instrument gates with provenance, without biological interpretation."""
import csv
import hashlib
import json
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import provenance
from elai_amy_gate import OUT


def main():
    metric = json.loads((OUT / 'gate_metrics.json').read_text())
    panel = json.loads((OUT / 'panel.json').read_text())
    runs = sorted((json.loads(p.read_text()) for p in OUT.glob('*.run.json')), key=lambda r: r['prefix'])
    with (OUT / 'samples.tsv').open() as f:
        samples = list(csv.DictReader(f, delimiter='\t'))
    groups = [s['group'] for s in samples]
    counts = {g: groups.count(g) for g in set(groups)}
    g = metric['gates']; a = g['known_structure']; s = g['stability']; v = g['virtual_sanity']
    hashes = json.loads((OUT / 'repeat_before_hashes.json').read_text())
    repeat = {p: hashlib.sha256(Path(p).read_bytes()).hexdigest() == h for p, h in hashes.items()}
    (OUT / 'repeat_validation.json').write_text(json.dumps(dict(files=len(repeat), identical=sum(repeat.values()), results=repeat), indent=2)+'\n')
    sizes = {p.name: p.stat().st_size for p in OUT.iterdir() if p.is_file()}
    (OUT / 'artifact_sizes.json').write_text(json.dumps(dict(total_bytes=sum(sizes.values()), max_file_bytes=max(sizes.values()), files=sizes), indent=2)+'\n')
    text = f'''# ELAI AMY fit results — gated STOP

[observation] **Verdict: {metric['verdict']}. Known-structure and stability gates FAILED.** Coverage, virtual-site arithmetic/integration, and exclusion gates passed. No allele-frequency analysis, temporal analysis, biological state interpretation or claim about AMY was performed. This is a bounded pilot on real imputed genotypes, not a full-panel validation.

## Input and instrument

[observation] Immutable source dimensions: {panel['source_dimensions'][0]:,} variants × {panel['source_dimensions'][1]:,} samples. Streamed {panel['pvar_records_scanned']:,} records, stopping at chr3. Exact AADR-1240K annotation/.ind matches, BP>0, deduplicated by annotation Individual ID: {panel['eligible_dated_unique']:,} eligible; deterministic random cohort {panel['fitted_n']} dated ancients, no modern references. {len(counts)} published groups ({sum(c == 1 for c in counts.values())} singleton groups); {len(set(s['country'] for s in samples))} country labels, {sum(s['country'] in ('..', '') for s in samples)} unknown countries. These sparse labels can inflate unadjusted NMI; NMI alone is not evidence of recovery.

[observation] Chr1 95,000,000–115,000,000; inclusive core 103,915,000–104,327,000 excluded. Fixed 2-kb candidate thinning gave {panel['candidates']:,} biallelic SNP candidates; cohort hardcall MAF≥0.05, call rate≥0.95 gave **{panel['markers']:,} fitted markers**, not the aspirational tens of thousands. Minimum retained MAF {panel['min_marker_maf']:.8f}, minimum call rate {panel['min_marker_call_rate']:.6f}. Alternate-marker thinning {panel['thin_markers']} markers. Thinning before QC often selected rare variants; we did not revise the preregistered rule after seeing results. Independent chr2 PC panel {panel['pc_markers']} SNPs (200-kb candidate spacing), first five standardized-genotype SVD PCs. This sparse PC panel and bounded cohort are limitations.

[observation] BIMBAM unphased diploid headers N/M/IND, allele strings mapped from native ALT hardcalls (not rounded fractional dosages). The imputed data's complete calls do not establish original aDNA coverage or imputation accuracy. Input panel.inp: {sizes['panel.inp']:,} bytes; panel.pos: {sizes['panel.pos']:,}; thin.inp: {sizes['thin.inp']:,}; thin.pos: {sizes['thin.pos']:,}. No matrix persisted over 500 MB. Full primary --pkk streamed {g['coverage']['full_export_rows']:,} rows through a FIFO; summary checks verified every numeric sample index, fitted marker ID and ordered row sequence. Upstream leaves individual_id blank despite the manual's optional IND header; samples.tsv plus unchanged input-column order is the authoritative index-to-sample mapping. Only virtual-site tensors, small composition summaries, parameters and logs retained.

[observation] Patched binary /home/erikg/elai-fork/src/elai, no rebuild; operator identifies its public fork as github.com/ekg/ELAI, default branch lower-state-posterior. Actual binary SHA256 is recorded in provenance, without relying on an unverified branch identity. Manual-derived settings: -p 1 (cohort), -C 1 (requested linked limit; upper state is constant), -c 5 (manual 5*C), sensitivities 3/8, -w 0 -s 20 (manual recommended exact EM 20–50), -mg 10 fixed prior scale, -R recorded seeds. pk1 COPY COUNTS divide by 2; ordered pkk sums to 1; q sums to 1 only on k1≤k2. The final likelihood still improves after 20 steps (primary last increment 119.681); convergence was not established. We do not tune EM or retry toward acceptance after observing gate failure.

## Numeric gates

| Gate | Observed | Preregistered requirement | Result |
|---|---|---|---|
| (a) Known structure | composition silhouette {a['silhouette']:.6f}; independent PC silhouette **{a['pc_silhouette']:.6f}**; group NMI {a['group_nmi']:.6f}; country NMI {a['country_nmi']:.6f}; {a['nonempty_assignments']} nonempty assignments | silhouette and PC silhouette >0.10; group NMI >0.10; country NMI >0.05 | **FAIL: PC silhouette** |
| (b) Stability | 8/8 comparisons completed; see below | seed/thinning ARI≥0.80, cross-K ARI≥0.60, same-K L1≤0.20, virtual agreement≥0.80, cross-K normalized Gram difference≤0.30 | **FAIL** |
| (c) Coverage | {g['coverage']['usable_individuals']}/{g['coverage']['fitted_n']} usable (100%); {g['coverage']['dated_virtual_output']} dated ancients in complete site output; {g['coverage']['full_export_rows']:,} complete all-marker rows | ≥95% individuals with ≥95% calls; ≥100 dated exports; all-marker tensor complete | PASS |
| (d) Virtual sanity | max ordered sum error {v['ordered_sum_error']:.3g}; q symmetry error {v['q_symmetry_error']:.3g}; triangular q error {v['q_triangle_sum_error']:.3g}; linked/unlinked max change **{v['bridge_max_change']:.9f}** | ≤1e-10, ≤1e-12, ≤1e-10 respectively; finite probabilities [0,1]; change >1e-6 | PASS |
| (e) Exclusion | 0 builder/retained core markers, 0 PC overlap; site 104121000 inside core; left 103783823 < site < right 104341989 | zero excluded markers; inside core; same-chromosome flanks | PASS |

[observation] Eight source-named helper regressions pass (tests.log); independent artifact validation confirms all nine executions, registration before fitting, exact zero-core-marker exclusion, complete 128-individual exports, q formula error 0, and exactly one edited control token (independent_validation.json).\n\n[observation] Other structure diagnostics (no minimum date correlation preregistered): date-bin NMI {a['date_bin_nmi']:.6f}, composition-predicted BP R² {a['date_composition_r2']:.6f}; published-group/country PC silhouettes {a['group_pc_silhouette']:.6f}/{a['country_pc_silhouette']:.6f}, composition-space silhouettes {a['group_composition_silhouette']:.6f}/{a['country_composition_silhouette']:.6f}. These do not support an independent-PC structure pass. Earlier painter silhouette<0.10/group NMI<0.10/country NMI<0.05 was a failed precedent; beating its NMI values is insufficient here.

### Stability numbers

| Comparison | ARI | Matched mean L1 | Virtual argmax agreement | Result |
|---|---:|---:|---:|---|
'''
    for r in s['seed']:
        text += f"| K={r['k']}, seeds 20261010/20261011 | {r['ari']:.6f} | {r['mean_l1']:.6f} | {r['virtual_argmax_agreement']:.6f} | FAIL |\n"
    for r in s['thinning']:
        text += f"| K=5, original/alternate thinning | {r['ari']:.6f} | {r['mean_l1']:.6f} | {r['virtual_argmax_agreement']:.6f} | FAIL |\n"
    text += '\n| Cross-K (seed) | ARI | Normalized Gram difference | Result |\n|---|---:|---:|---|\n'
    for r in s['cross_k']:
        text += f"| K={r['k_pair'][0]}/{r['k_pair'][1]} ({r['seed']}) | {r['ari']:.6f} | {r['gram_difference']:.6f} | FAIL |\n"
    text += f'''
[observation] Composition-based Hungarian matching is applied only at the same K. Cross-K comparison uses ARI and normalized sample co-composition Gram matrices, not arbitrary label correspondence. Low Gram differences with poor ARI show why marginal similarity alone cannot pass stability. Virtual agreement is the lower marginal argmax, not an inferred phased haplotype. Prior painter's failed 43% target agreement is not rescued: this pilot gives 7.8125% seed agreement at K=5 and 14.84375% thinning agreement.

## Exclusion and virtual control demonstrated

[observation] Builder and runner assert inclusive exclusion; the runner rechecks every actual retained snpinfo record in all nine successful executions, and checks the virtual locus/coordinate and complete K*K rows for all individuals. Virtual ordered-copy symmetry max error {v['ordered_symmetry_error']:.3g}. Frozen control preserves all EM-file lines and changes only rk at the right core-gap flank to 1 (fully unlinked); original fitted file is replayed with s=0 as the linked rounded-parameter baseline. The two outputs differ by {v['bridge_max_change']:.9f} maximum absolute ordered probability on this input, confirming linkage sensitivity, not biological accuracy of the missing-emission bridge. Parameters, edited control, commands and outputs are retained.

## Exact reproducible executions

[observation] Run from repository root with .venv/bin/python, OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1:
```
.venv/bin/python analysis/elai_amy_build.py
.venv/bin/python analysis/elai_amy_gate.py
.venv/bin/python analysis/elai_amy_gate.py --metrics-only
PYTHONPATH=analysis .venv/bin/python -m unittest analysis/test_elai_amy.py -v
.venv/bin/python analysis/elai_amy_report.py
```
The runner creates/deletes a scratch working directory because ELAI always writes relative output/. The exact binary commands below use persistent input paths; prefix uniquely identifies K/panel/seed. Successful upstream exit code is 1; finished sentinel, output completeness and numeric checks are also required.

'''
    for r in runs:
        text += f"[observation] {r['prefix']}: {r['markers']} retained markers, {r['elapsed_seconds']:.3f} seconds, exit {r['returncode']}.\n```sh\n{shlex.join(r['command'])}\n```\n\n"
    text += f'''## Execution audit, resources, validation and decision

[observation] Registration hash {panel['prereg_sha256']} is unchanged across runs; its filesystem mtime precedes every recorded final fit start. Three execution-only repairs are preserved: first parser assumed the manual's six-column snpinfo description, but this binary actually writes rs,pos,chr,major,minor,maf,miss,n0,n1,n2; second control failed because parameter files are line-oriented and flattening tokens breaks -rem; a repeat-audit assertion initially required textual IDs even though README explicitly warns BIMBAM IDs may be empty. Corrected parsing/line preservation and authoritative numeric-index auditing, with **no changes to cohort, markers, model settings or gates**. Historical failure logs/metrics are explicitly named first/control/pre_control. A complete final rerun reproduced {sum(repeat.values())}/{len(repeat)} site/composition file SHA256 hashes byte-for-byte.

[observation] Starting, pre-write and final df -h / show 71G free (rounded 100% use on shared 14T volume); final metric-time free {metric['disk_free_bytes']:,} bytes. Persistent lane output approximately {sum(sizes.values())/1e6:.2f} MB (max file {max(sizes.values())/1e6:.2f} MB), far below 5GB. Scratch removed and deletion recorded in scratch_cleanup.json; no raw inputs modified, no installs, no rebuild, no git mutation/staging. Final permitted read-only `git diff --cached --name-only` returned no entries (no staged files); no git command ran in elai-fork. Provenance uses scripts/provenance.py with git provider disabled before construction; small files and binary hashed, huge immutable pgen/pvar recorded by size/mtime rather than claiming fresh content hashes.

[inference] **STOP is the outcome, not a biological negative result.** Required before a new attempt: prospectively registered denser QC-then-thin markers, more representative/larger cohort and independent genome-wide PC validation, adequate EM convergence, and then renewed seed/K/thinning validation with the same fail-closed discipline. This lane does not establish that ELAI generally fails, nor that any AMY history is absent. No local-state biological summary or effective independent-state estimate is justified because gates failed. Temporal and allele-frequency analyses remain blocked.
'''
    result_doc = ROOT / 'docs/elai-amy-fit-results.md'
    result_doc.write_text(text)
    provenance.get_git_info = lambda: dict(commit='not queried', branch='not queried', dirty='not queried; operator prohibits git probes')
    with provenance.RunRecorder('elai-amy-gated-stop') as recorder:
        recorder.set_seed(20261010)
        for p in [ROOT/'docs/elai-amy-prereg.md', ROOT/'analysis/load_akbari_geno.py', ROOT/'scripts/provenance.py', *ROOT.glob('analysis/*elai_amy*.py'), Path('/home/erikg/elai-fork/src/elai'),
                  Path('/home/erikg/elai-fork/README-PATCH.md'), Path('/home/erikg/elai-fork/elai-manual.tex'), Path('/home/erikg/elai-fork/elai-manual.pdf'),
                  *ROOT.glob('data/raw/aadr-v66.p1/v66.p1_1240K*.anno'), *ROOT.glob('data/raw/aadr-v66.p1/v66.p1_1240K*.ind')]:
            recorder.record_input(p)
        for p in ROOT.glob('data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.*'):
            if p.stat().st_size < 500_000_000:
                recorder.record_input(p)
            else:
                recorder.manifest['inputs'][str(p)] = dict(bytes=p.stat().st_size, mtime_ns=p.stat().st_mtime_ns, hash='not computed: bounded I/O')
        recorder.record_parameter('exact_elai_commands', [r['command'] for r in runs])
        recorder.record_parameter('prereg_mtime', (ROOT/'docs/elai-amy-prereg.md').stat().st_mtime)
        recorder.record_parameter('gates', metric)
        for p in OUT.iterdir():
            if p.is_file():
                recorder.record_output(p)
        recorder.record_output(result_doc)
    print(result_doc)


if __name__ == '__main__':
    main()
