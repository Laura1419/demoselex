# Research Journal

This journal records the progress and decisions of the trajectory-replication research project.

**Newest entries first. Format:**
```
## YYYY-MM-DD HH:MMZ — <task-id> — <one-line>
- what was done
- what was found
- blockers / next
```

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