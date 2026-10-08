#!/usr/bin/env python3
"""Create the living-document skeleton (docs sections A-J + support docs).

Idempotent: existing files are left untouched unless --force.
"""
from __future__ import annotations

import os
import sys

ROOT = "/home/erikg/demoselex"

DOCS: dict[str, str] = {
    "docs/A-scientific-summary.md": """# A. Scientific summary

> Owner task: recon-akbari / synth-assessment. Status: TODO.

Plain-language statement of the question and why it matters. Write for a
quantitative geneticist who does not know this project.

## The two scenarios
## Why the population-structure correction may hurt
## What would change our minds
""",
    "docs/B-literature-map.md": """# B. Literature map

> Owner task: recon-lit. Status: TODO.

| Method | Unit of analysis | Structure representation | Nuisance | Evidence of selection | Per-population selection? | Repeated direction across histories? | Adaptable here? |
| --- | --- | --- | --- | --- | --- | --- | --- |
| (fill) | | | | | | | |

## Novelty statement
""",
    "docs/C-akbari-reconstruction.md": """# C. Akbari et al. 2026 reconstruction

> Owner task: recon-akbari. Status: TODO.

## Paper and identifiers
- Akbari et al. 2026, Nature 654:419-428, doi:10.1038/s41586-026-10358-1
- Preprint: bioRxiv 2024.09.14.613021 / PMC11419161
- ENA: PRJEB106907

## Sample and scope
## Pipeline
## The model (likelihood written out)
## Selection coefficient: definition and units
## Locus-specific variance
## Calibration and filters
## What is downweighted under time-GRM collinearity
## Released data
""",
    "docs/D-method.md": """# D. Proposed method

> Owner task: lean-* + sim-stat. Status: TODO.

## Definitions
## Stage A: trajectory inference
## Stage B: replicated directional change
## Statistics (Model 0 .. Model 3)
## Null calibration
## One-paragraph plain-language summary
> "We inferred major population trajectories from genome-wide variation independently
> of each tested locus. We then estimated allele-frequency change through time along
> each trajectory and tested whether the direction and magnitude of change were
> concordant across genetically distinct trajectories, while accounting for their
> non-independence. Significance was calibrated against matched genomic control loci
> and demographic simulations."
""",
    "docs/E-data-inventory.md": """# E. Data inventory

> Owner task: recon-data. Status: TODO.

| Source | Accession/URL | Format | Size | License | Access | Download command | Fallback |
| --- | --- | --- | --- | --- | --- | --- | --- |
""",
    "docs/F-reproducible-code.md": """# F. Reproducible code

> Owner tasks: all. Status: TODO.

Inventory of scripts, what each consumes/produces, and the exact command to run it.
""",
    "docs/G-validation-report.md": """# G. Validation report

> Owner tasks: sim-sweeps, locus-canonical, traj-validate. Status: TODO.

## Decision 1 — trajectory stability
## Decision 2 — null calibration
## Decision 3 — power vs Akbari, FPR under migration
## Decision 4 — AMY and canonical loci
""",
    "docs/H-amy-analysis.md": """# H. AMY analysis

> Owner task: locus-amy. Status: TODO.

1. Does AMY copy number rise along multiple independent trajectories?
2. Are rises confined to replacements?
3. Do several structural haplotypes rise?
4. Do independent mutational backgrounds converge?
5. Does the evidence survive excluding close trajectories?
""",
    "docs/I-akbari-reanalysis.md": """# I. Akbari reanalysis

> Owner task: locus-akbari. Status: TODO.

## Four-quadrant classification
| Quadrant | Akbari | Rep  | Interpretation |
| --- | --- | --- | --- |
| Q1 | high | high | strongest consensus |
| Q2 | high | low  | single broad/demographic event |
| Q3 | low  | high | **candidate over-correction** |
| Q4 | low  | low  | weak/null |

## Low-Akbari / high-replication loci in depth
""",
    "docs/J-final-assessment.md": """# J. Final assessment

> Owner task: synth-assessment. Status: TODO.

Verdict (choose one, with evidence chain):
- [ ] strongly supported
- [ ] partly supported
- [ ] unsupported
- [ ] dependent on particular forms of selection
""",
    "docs/DECISIONS.md": """# Decision log (plan §20)

| # | Decision | Status | Evidence | Date |
| --- | --- | --- | --- | --- |
| 1 | Can stable ancestry trajectories be inferred? | OPEN | | |
| 2 | Is the replication statistic calibrated under neutral demography? | OPEN | | |
| 3 | Does the method beat Akbari under simulated repeated selection without excess FPR? | OPEN | | |
| 4 | Do AMY and canonical loci behave sensibly? | OPEN | | |
| 5 | Is there a meaningful low-Akbari/high-replication class? | OPEN | | |
""",
    "docs/PROVENANCE.md": """# Provenance policy

Every run records: git commit, dirty flag, command line, software versions, input
file hashes, parameters, random seed, outputs. Written to
`analysis/runs/<UTC timestamp>-<slug>/manifest.json` by `scripts/provenance.py`.

Raw data under `data/raw/` is immutable. All derived data is regenerable.
""",
    "docs/JOURNAL.md": """# Research journal

Append-only. Newest entries first. Format:

```
## YYYY-MM-DD HH:MMZ — <task-id> — <one-line>
- what was done
- what was found
- blockers / next
```
""",
}

HEADERS = {
    "docs/A-scientific-summary.md": None,
}


def main() -> int:
    force = "--force" in sys.argv
    for rel, body in DOCS.items():
        path = os.path.join(ROOT, rel)
        if os.path.exists(path) and not force:
            print(f"keep: {rel}")
            continue
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(body)
        print(f"write: {rel}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
