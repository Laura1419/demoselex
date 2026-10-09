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
