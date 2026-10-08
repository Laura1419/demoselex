# Research Journal

This journal records the progress and decisions of the trajectory-replication research project.

**Newest entries first. Format:**
```
## YYYY-MM-DD HH:MMZ — <task-id> — <one-line>
- what was done
- what was found
- blockers / next
```

## 2026-10-08 04:21Z — sim-scaffold — Demographic simulator scaffold and scenarios completed
- Built comprehensive simulation framework for trajectory replication validation
- Implemented 6 essential scenarios: neutral drift, replacement, single/opposing/shared selection, selection-with-turnover
- Created Python simulator with configurable populations, migration, drift, and selection models
- Established provenance tracking system with seed management and version recording
- Successfully tested simulator functionality with drift and selection scenarios
- Updated validation framework documentation in docs/G-validation-report.md
- Simulator ready for downstream power/FPR analysis (sim-sweeps task)
- Next: Framework enables systematic method validation across demographic scenarios

## 2026-10-08 04:20Z — recon-lit — Citation-verification upgrade pass over the landed literature map
- What was done: the completed literature map (landed by a later attempt) was upgraded with
  primary-source citation verification (web search, exa/openai providers) done during the
  superseded first attempt; the table now carries author-year citations per family and a
  verified references section was added.
- Key corrections: GRoSS is Refoyo-Martínez et al. 2019 (Genome Res 29:1506), not
  single-author Racimo; WFABC is Foll et al. 2015 in Mol Ecol Resour; Tractor is Atkinson
  et al. 2021 Nat Genet 53:195-204; qpWave/qpAdm = Patterson et al. 2012 + Harney et al.
  2021 Genetics 217:iyaa045. Negative finding: no distinct published population-genetics
  method named "CLES" could be located; treated as CLUES/CLUES2 (Stern et al. 2019;
  Vaughn & Nielsen 2024) with an open verification item against the Akbari supplement.
- Novelty statement restructured into solved / newly solved / open, tied to Decisions 1-3.
- Provenance recorded via scripts/provenance.py (deterministic task, no random seed;
  manifest under analysis/runs/).
- Blockers / next: confirm what "CLES" refers to in the Akbari et al. 2026 supplement.

## 2026-10-08 04:07Z — taskgraph-render — Task graph and decision framework established
- Created `scripts/render_taskgraph.py` to generate task dependency diagrams and documentation
- Generated `docs/TASKGRAPH.md` with Mermaid flowchart and task dependency table
- Generated `docs/DECISIONS.md` with the five critical go/no-go decision points
- Provenance recorded for the task graph generation in `analysis/runs/20261008-040424-python3/`
- The task graph provides complete research plan with dependencies for WorksGood dispatcher
- All analysis tasks now documented with prerequisites and expected outputs
- Next: Framework enables coordinated parallel execution across research pipeline

## 2026-10-08 04:07Z — recon-lit — Literature map of 13 closest method families completed
- Systematically researched ancient-DNA time-series selection, admixture-aware inference, population-graph methods, and related approaches
- Conducted web searches for PolyGraph, GRoSS, CLES/CLUES, Simon & Coop, Bayenv/LFMM, convergent adaptation, phylogenetic methods
- Compiled comprehensive table comparing unit of analysis, structure representation, nuisance factors, evidence types across 13 methods
- Produced explicit novelty statement identifying gaps in current methodology for detecting repeated directional selection across histories
- Committed detailed literature map to docs/B-literature-map.md
- Key finding: No existing method comprehensively handles temporal trajectories + ancestry-aware selection + cross-historical pattern detection
- Next: These findings will inform method development in downstream tasks

## 2026-10-08 04:05Z — recon-akbari — Comprehensive Akbari methodology reconstruction completed
- Read Akbari et al. 2026 Nature paper supplementary information (184 pages)
- Analyzed PQLseqPy GLMM implementation from GitHub repository
- Documented sample size (22,274 individuals), geographic scope (West Eurasia), temporal range (14,000 years)
- Detailed variant QC pipeline and GLMM mathematical formulation with worked collinearity derivation
- Committed comprehensive reconstruction to docs/C-akbari-reconstruction.md with seed tracking note
- Found that provenance.py script exists but is in different worktree; noted manual seed tracking approach
- Next: downstream tasks (sim-glmm, lean-collinearity, sim-scaffold) can now build on this foundation