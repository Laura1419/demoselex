# Provenance and Reproducibility

This document describes the provenance and reproducibility framework for the trajectory-replication project.

## Overview

Every analysis run records a comprehensive manifest containing:
- **Git commit** - exact code version
- **Command line** - execution parameters
- **Software versions** - dependencies and system info
- **Input file hashes** - content verification
- **Parameters** - configuration and random seed
- **Outputs** - produced files and their hashes

## Manifest Format

Manifests are written to `analysis/runs/<timestamp>-<slug>/manifest.json` with the following structure:

```json
{
  "run_id": "20241008-1430-trajectory-analysis",
  "timestamp": "2024-10-08T14:30:00Z",
  "git": {
    "commit": "a1b2c3d4e5f6...",
    "branch": "main",
    "dirty": false
  },
  "command": ["python", "analysis/trajectory.py", "--config", "config.json"],
  "versions": {
    "python": "3.11.0",
    "pandas": "2.1.0",
    "numpy": "1.24.0"
  },
  "inputs": {
    "data/genotypes.vcf": "sha256:abcd...",
    "config.json": "sha256:efgh..."
  },
  "parameters": {
    "seed": 42,
    "iterations": 1000,
    "threshold": 0.05
  },
  "outputs": {
    "trajectories.csv": "sha256:ijkl...",
    "plot.png": "sha256:mnop..."
  }
}
```

## Usage

### Library Interface

```python
from scripts.provenance import RunRecorder

# Start recording a run
with RunRecorder("trajectory-analysis") as recorder:
    # Record parameters
    recorder.record_parameter("seed", 42)
    recorder.record_parameter("iterations", 1000)

    # Record input files
    recorder.record_input("data/genotypes.vcf")
    recorder.record_input("config.json")

    # Execute your analysis
    result = run_analysis()

    # Record outputs
    recorder.record_output("trajectories.csv")
    recorder.record_output("plot.png")
```

### CLI Tool

```bash
# Record a one-off command
python3 scripts/provenance.py record \
    --slug "quick-test" \
    --input data/genotypes.vcf \
    --input config.json \
    --parameter seed=42 \
    --command "python analysis/trajectory.py --config config.json"

# Check reproducibility of a run
python3 scripts/provenance.py verify analysis/runs/20241008-1430-trajectory-analysis/manifest.json
```

## Principles

1. **Raw data is immutable** - Never modify files in `data/raw/`
2. **Derived data is regenerable** - All derived data can be recreated from raw + code
3. **Seed all randomness** - Every random process gets a deterministic seed
4. **Hash all inputs/outputs** - Track content, not just filenames
5. **One manifest per run** - Each execution creates its own record

## Integration with Pipeline

The provenance harness is integrated into all analysis pipeline steps:

```python
# Example pipeline step
from scripts.provenance import run_with_provenance

def analyze_trajectories():
    return run_with_provenance(
        slug="trajectory-analysis",
        inputs=["data/genotypes.vcf", "config.json"],
        parameters={"seed": 42, "iterations": 1000},
        command_fn=lambda: real_analysis_function()
    )
```

## Validation

To validate the system:

1. Run the demo: `python3 scripts/provenance.py demo`
2. Check the generated manifest contains expected fields
3. Verify git commit, versions, and file hashes are recorded
4. Test reproducibility by verifying a recorded run

The system is considered working when a demo run produces a complete manifest with all required metadata.

## Random Seed Management

All random processes must use seeded random generators:

```python
import random
import numpy as np

seed = recorder.get_seed()  # Returns deterministic seed from run context
random.seed(seed)
np.random.seed(seed)
```

This ensures that rerunning with the same seed produces identical results, allowing validation of reproducibility.

## Third-pass prioritisation manifests — 2026-10-09 (append only)

[observation] Completed manifests: analysis/runs/20261008-2359-prioritization/manifest.json and analysis/runs/20261009-0002-prioritization/manifest.json. The latter includes exact pvar-cache index verification for10,642 panel variants. Generated with scripts/provenance.py RunRecorder, seed20261009; git probes monkeypatched OFF to obey the operator's no-git requirement. Small source/input/output artifacts are SHA256 hashed; immutable large genotype/cache inputs are recorded by path and byte size rather than rehashed. The first attempt stopped before ranking; its failure and repair are preserved in logs/preregistration. No pristine external calibration is claimed.

## Appended 2026-10-09 — block-survey provenance

[observation] Annotation source: `https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz`, retrieved 2026-10-09T02:19:41.016138+00:00, verified gzip CRC; compressed SHA256 `28e7fe22cf1923a6887c7b60e6004c75987f810b8374d4d396021a34ac0410d1`. Immutable download and metadata reside in `data/raw/annotation/`. Coordinates are zero-based half-open transcript intervals in GRCh37/hg19; source symbols and strand are retained.

[observation] `analysis/block_survey.py` uses `scripts.provenance.RunRecorder`, overriding only `get_git_info` to record “operator prohibits all git commands” before constructing any recorder. It hashes all four shard TSVs, pre-registration, code, annotation and metadata, PDF, workbook, QC cache, saved simulation/metrics inputs, and derived survey outputs. Initial bounding artifacts are retained under `data/derived/block_survey_initial_bounding*` to audit the approved interval-union correction. The results doc is append-only; the latest appended survey supersedes historical annotation attempts. Run/validation manifests are in `analysis/runs/*block-survey*/manifest.json`.

## Two-channel provenance — appended 2026-10-09T02:59Z

[observation] Registration: docs/two-channel-prereg.md persisted 02:49:12Z before first result approximately 02:51Z. RunRecorder from scripts/provenance.py records seed 20261009, software versions, commands, ancestry-relatedness input, source hashes and output hashes. All git probes are disabled BEFORE constructing the recorder, per operator rule. AMY bootstrap uses seed20261010; each correlated drift scenario uses4096 draws. Full SNP/AMY/audit modes were each repeated; final production manifests are analysis/runs/20261009-0256-two-channel-genomewide/manifest.json, analysis/runs/20261009-0254-two-channel-amy/manifest.json and analysis/runs/20261009-0254-two-channel-audit/manifest.json. Earlier manifests remain intact.

[observation] SNP source panel/priority/sample-metadata TSVs are hashed for all eight directories; cached genotype arrays and the674,584,607-byte immutable statistics release are recorded by path/byte size rather than rehashed. AMY raw CN and archive are SHA256 hashed, as are all generated channel tables. This bounded-I/O exception is explicit, not a content-verification claim. The full two_channel.tsv hash in the final production manifest is e5bcc37470381a4917131f7220cbf5f6b77fc00022b688412e5fa5edd743d36c. Raw inputs were not edited. No git command/staging/commit was run.

[observation] Validation provenance additionally records five new tests/seven existing tests, exact artifact coverage, repeated-run identical hashes for amy_summary.json/all eight SNP tables, source/docs and generated block/validation outputs. The final validation manifest under analysis/runs/*two-channel-validation-final/ captures final logs/docs; the earlier validation record remains intact (its test log was subsequently re-run, changing elapsed time only). No calibrated LD-honest inference is inferred from these software checks.

## CN-as-SNP provenance — appended 2026-10-09

[observation] Preregistration was created before biological fitting; its pre-fit file mtime is preserved in the first-run summary and independently checked against run start. A clearly marked secondary numerical-diagnostic amendment was appended after primary fits; it is not presented as preregistered before those fits. Primary runs use scripts/provenance.py RunRecorder under analysis/runs/*cn-as-snp*/; all git metadata probes are disabled before constructing the recorder, per the task's no-git rule. The first summary/statistics are retained as first_run_* snapshots; after correcting only a misleading completion-status label, the second primary run exactly reproduced the biological result/null rows.

[observation] Secondary diagnostic and unrestricted-null runs plus the final validation/report session have dedicated manifests. Source, registration, archive, sample metadata, Table S1, outputs and logs are SHA256-recorded as appropriate. Large immutable pgen/pvar/statistics inputs are recorded by byte size/mtime rather than newly asserted content hashes. No git identity or clean-index assertion is fabricated; no git commands/staging were performed. All six saved LOCO matrices were independently recomputed from genome-wide sample-matched genotypes. Seeds=20261009; 199 attempted permutations per available locus/model, with numerical failures explicit.

## APPENDED full-panel haplotype-background provenance — 2026-10-09

[observation] Code `analysis/haplotype_background.py` and append-only reporter use `scripts/provenance.py` RunRecorder with git metadata probes disabled BEFORE construction. Initial-method manifests: `analysis/runs/20261009-1737-haplotype-background/manifest.json`, `analysis/runs/20261009-1740-haplotype-background/manifest.json`; corrected manifests: `analysis/runs/20261009-1749-haplotype-background/manifest.json`, `analysis/runs/20261009-1751-haplotype-background/manifest.json`. Independent baseline/scan checks and the report have separate17:43/17:50/17:52 manifests. Initial outputs, including original coordinates and assignments, are retained under `data/derived/haplotype_background/invalid_exclusion_run/`. The original preregistration bytes retained as `prereg.original.md` match original SHA256 ca49cf2edc71b4d85554a49cb8545c921766d167e62bafdf1063385b27f9772e; its original pre-fit mtime is in first-run records. Post-results annotation correction is explicitly appended, with no K/gate tuning.

[observation] All16,339 full-panel rows participate at every locus. Raw psam/annotation/archive/source/preregistration/small artifacts are hashed; large immutable pgen/pvar/released statistics are recorded by size/mtime rather than asserting new content hashes. Corrected repeat statistics are identical apart from timing/RSS. Final validation manifest additionally captures current source/test/doc hashes and stable completed logs; earlier run manifests are historical snapshots, and files that changed under correction/reporting must not be compared to earlier snapshot hashes as if unchanged. No raw writes and no git command, staging or commit were performed. The index was not inspected under the no-git rule.

## APPENDED Akbari peer-review evaluation provenance — 2026-10-09T20:02:50+00:00

[observation] Source data/raw/akbari/41586_2026_10358_MOESM4_ESM.pdf SHA2568c94cddd5cd6c851403f68a4c3168dedfe5f5d29d1a973a808c0a0960365d108 verified;61 pages extracted with .venv/bin/python and pypdf. Derived physical-page JSON/text,264-entry quote ledger and validation summary are under data/derived/akbari_peer_review/. Three analysis scripts preserve extraction, one-shot append generation and repeatable source/section/table/absence validation. A scripts/provenance.py RunRecorder manifest records source/scripts/final report/derived hashes; git provider disabled BEFORE construction under the operator’s no-git rule. No raw writes, git command or staging; existing repository index state not inspected. Figures are not numerically reconstructed, and rebuttal numerical claims are not independently adopted as results.

## APPENDED ELAI AMY instrument-gate provenance — 2026-10-10

[observation] analysis/elai_amy_build.py streams pvar/pgen into the preregistered bounded cohort/panel; analysis/elai_amy_gate.py records every exact binary invocation, seed, start epoch, preregistration mtime, retained markers, arithmetic checks and FIFO coverage. Registration SHA256 ae4ae3e2dd0942abd40e43bbcfbd3b618e046803d0d2dde48bf77790c614d460 predates all final fits and is unchanged. Final runs: seven fits (K3/5/8, seeds20261010/20261011; alternate thinning K5) plus two s=0 controls. Full primary tensor is streamed, never persisted; only small site exports/composition summaries/parameters/logs retained. The control edits exactly token2872 (gap rk=0.3746716→1), preserving upstream line layout. All16 retained site/composition files reproduce byte-identically across complete reruns.

[observation] scripts/provenance.py RunRecorder is used by analysis/elai_amy_report.py with git provider disabled before construction. Manifests under analysis/runs/*elai-amy-gated-stop*/ record binary/manual/patch README/source/psam/AADR hashes, exact commands and gate metrics; pgen/pvar use immutable size/mtime metadata rather than newly asserted full hashes. Example completion snapshot: analysis/runs/20261010-0145-elai-amy-gated-stop/manifest.json; earlier snapshots preserve execution/audit history and are not claims that subsequently regenerated report/log bytes remain identical. Independent_validation.json verifies registration order, exclusion, dated count, unordered formula and single-token control edit; tests.log records eight passing regressions. Three IO/assertion-only repairs are disclosed in the results, without changing the preregistered design or thresholds.

[inference] Scientific verdict STOP: known-structure independent-PC silhouette=-0.057757 and unstable seeds/K/thinning; coverage, exclusion and virtual arithmetic/integration pass. No AMY biological interpretation or allele-frequency result exists. Scratch is removed, derived outputs ~5.4MB, volume free space71G. No raw edits, ELAI repo writes, installs, rebuilds or git mutation/staging/commit; permitted final read-only git diff checks showed no staged files.
