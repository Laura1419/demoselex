# MISSION — master specification

This file is the authoritative statement of the objective. Work items in `.wg`
implement it. If a task and this file conflict, this file wins until deliberately
amended (amendment = commit with `MISSION:` prefix).

## Objective

Determine whether *replicated directional allele-frequency change across independently
inferred population trajectories* adds real information beyond the population-structure
correction used by Akbari et al. (2026), and characterize the set of loci where it does.

## The idea in one paragraph

Ancient-DNA selection scans correct for population structure because migration and
turnover produce allele-frequency changes that mimic selection. But that correction may
also remove real selection when selection is correlated with population history. Instead
of treating structure purely as nuisance, infer a small set of distinct temporal
population trajectories from genome-wide data, then ask whether a candidate allele
repeatedly changes in the same direction *within* multiple trajectories. A locus rising
independently along several genetically distinct trajectories is stronger evidence than a
locus whose entire change comes from one high-frequency population replacing a
low-frequency one.

## Required outputs (§17)

Maintain a living research document with sections:

- **A. Scientific summary** — plain language, why it matters.
- **B. Literature map** — table of related methods and how they differ.
- **C. Akbari reconstruction** — exact model + available data.
- **D. Proposed method** — equations + plain-language interpretation.
- **E. Data inventory** — files, repositories, metadata, access requirements.
- **F. Reproducible code** — download/QC/trajectory/locus/sim/figures/reanalysis.
- **G. Validation report** — simulations + known loci.
- **H. AMY analysis** — motivating case study.
- **I. Akbari reanalysis** — reported loci and, if feasible, genome-wide.
- **J. Final assessment** — strongly / partly / not supported / dependent on selection form.

## Decision points (§20) — explicit go/no-go

1. Can stable ancestry trajectories actually be inferred? If no → reformulate.
2. Under neutral demographic simulations, is the replication statistic calibrated?
   If no → fix before touching biological loci.
3. Does the method recover simulated repeated selection better than Akbari without
   greatly increasing false positives? If no → the method has not solved the problem.
4. Do AMY and canonical loci behave sensibly? If no → investigate before scaling.
5. Does the Akbari reanalysis reveal a meaningful low-Akbari/high-replication class?
   If no → conclude the hypothesised over-correction appears limited under this test.

Each decision point is a WG task whose completion is reported in `docs/DECISIONS.md`.

## Guardrails (§14, §18)

- Population expansion ≠ locus-specific selection. Require one or more of: within-trajectory
  change; repeated same-direction change across distinct trajectories; excess change beyond
  ancestry-mixture expectation; independent molecular mutations → same functional state;
  independent geographic/cultural replicates.
- Never claim over-correction without simulated or empirical controls showing a class of
  true signals that Akbari systematically suppresses.
- Quantify trajectory non-independence before treating trajectories as replicates.
- Do not force the project toward the original hypothesis. "Unsupported" is a valid result.

## Formalization requirement

Everything that can be stated mathematically is stated in Lean 4 (`lean/`). Definitions:
temporal trajectory, allele-frequency path, direction concordance, effective number of
independent trajectories. Theorems to attempt: (i) replacement/null ⇒ concordance
statistic has the exchangeable null distribution; (ii) time–GRM collinearity ⇒ loss of
identifiability of the GLMM slope (variance inflation); (iii) replication ⇒ detectable
power separation from replacement under stated conditions. Proofs may be partial; all
models must be *stated*, with honest `sorry`/`axiom` accounting in `lean/AXIOMS.md`.

## Continuous operation

The WorksGood task graph is the persistent plan. The dispatcher (`wg service start`)
runs ready tasks. Every task: log provenance, commit artifacts, update the relevant doc
section and `docs/JOURNAL.md`, and only then `wg done`.

## OBJECTIVE AMENDMENT: prioritisation — 2026-10-08

[observation] Operator directive for the third pass: “The scan is a PRIORITIZER for downstream phenotype work, not a discovery claim”; “False positives are EXPECTED AND ACCEPTABLE.” The requested output is a ranked ancient-DNA candidate list for modern association/phenotype follow-up.

[inference] For this purpose, ranking quality (positive-control enrichment, precision@N, recall, stability and drift-only follow-up cost) supersedes FDR gating. Demographically carried loci are down-weighted, not excluded. Earlier hypothesis-test framing remains valid for its own scientific question, but is NOT the current prioritisation objective. All earlier text and conclusions are retained.
