# Replicated Directional Change Across Population Trajectories

Testing whether population structure can be *used* as evolutionary replication rather
than only corrected away, and comparing the resulting statistic against the
Akbari et al. (2026) ancient-DNA selection scan for West Eurasia.

- Paper: Akbari, A. et al. *Ancient DNA reveals pervasive directional selection across
  West Eurasia.* Nature 654, 419–428 (2026). doi:10.1038/s41586-026-10358-1
- Preprint: bioRxiv 2024.09.14.613021 / PMC11419161
- Supplementary: `10.1038/s41586-026-10358-1` ESM (sections S1–S6)
- Data: ENA `PRJEB106907`; AADR / Reich lab releases; AMY (Bolognini et al. 2024)

## Core scientific question

When genetically different populations are followed through time, does the same locus
repeatedly move in the same direction?

- **Scenario A (replacement):** pooled frequency shifts because population B (f=0.7)
  replaces population A (f=0.1). Not by itself evidence of selection.
- **Scenario B (replicated directional change):** several genetically distinct
  trajectories each move in the same direction. Repeated movement along distinct
  histories is evidence that a single ancestry replacement cannot explain.

Primary hypothesis: *some loci whose signals are weakened by population-structure
correction nevertheless show repeated directional allele-frequency change across
multiple independently inferred population trajectories.* Tested, not assumed.

## Repository layout

| Path | Purpose |
| --- | --- |
| `docs/` | Living research document (§A–§J of the project plan) |
| `lean/` | Lean 4 + Mathlib formalization of the core objects and theorems |
| `sim/` | Forward/coalescent simulations, Akbari-GLMM reimplementation, power/FPR sweeps |
| `analysis/` | Data QC, trajectory inference, locus analysis, Akbari reanalysis |
| `scripts/` | Download / provenance / orchestration utilities |
| `fig/` | Regenerable figures |
| `data/raw/` | Immutable downloads (gitignored) |
| `data/derived/` | Reproducible intermediates (gitignored) |
| `.wg/` | WorksGood task graph (the tasking graph driving continuous work) |

## Modeling hierarchy (do not skip ahead)

- **Model 0** — trajectory-specific slopes of allele dosage ~ time.
- **Model 1** — weighted meta-analysis across trajectories.
- **Model 2** — covariance correction among trajectories (GLS / phylogenetic meta-analysis).
- **Model 3** — hierarchical shared-mean + heterogeneity model.

## Evidence standards

Distinguish **observation / inference / speculation**. Do not call a trajectory
independent without quantifying its genetic non-independence. Do not call a locus
selected merely because it changes monotonically. Do not claim Akbari "over-corrects"
without simulations or empirical controls. Do not treat archaeological labels as
ground truth.

## Reproducibility

Every data source, accession, download date, software version, parameter, filter,
seed, and exclusion is recorded in `docs/PROVENANCE.md` and per-run manifests under
`analysis/runs/`. Raw data immutable; all principal figures regenerable from code.

## Task graph

Work is decomposed into a dependency graph in `.wg` (WorksGood). See
`docs/TASKGRAPH.md` for the rendered graph and the decision points (§20 of the plan).

```bash
wg list --all          # list tasks
wg show <id>           # inspect a task
wg ready               # tasks with no unmet dependencies
wg service start       # run the dispatcher for continuous work
```
