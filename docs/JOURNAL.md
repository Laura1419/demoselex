# Research Journal

This journal records the progress and decisions of the trajectory-replication research project.

**Newest entries first. Format:**
```

## 2026-10-08 04:57Z — recon-data — Inventory handoff audit
- Rechecked the committed Section E access matrix against this task's validation: each requested source has an accession, exact URL, retrieval or status command, access status, and a fallback for controlled/embargoed records. No raw data or Lean files changed.
- The earlier provenance JSON recorded sources/seed/versions but was not in `scripts/provenance.py`'s verifiable format. Added a script-native manifest with seed 0 and hashed inventory/journal outputs; retained the historical manifests unchanged.

## 2026-10-08 04:39Z — recon-data — Primary-source file/access inventory
- Replaced the Section E placeholder with `docs/E-data-inventory.md`: exact Dataverse file IDs/byte sizes for Akbari statistics, imputed BCF/PLINK2 and AADR v66.p1 1240K/HO, ENA per-run manifest recipe, Akbari Nature supplements, and pinned AMY Zenodo/GitHub tables.
- Corrected previous attempt's unsupported assertions: the selection summary is public CC0, not subscription-only; ENA PRJEB106907 exposes aligned BAM/FASTQ rather than an asserted dosage VCF; AADR v66.p1 supersedes v66.0 and old Reich index link is unavailable.
- Access finding: Akbari's precise unpublished archaeology/coordinates XLSX is embargoed until 2041-01-01 (API 403), while study-used point dates and broad regions are public. Use public AADR `.anno` for overlapping samples only; never fabricate missing coordinates.
- Transfer caveat: public Dataverse API redirect followed by S3 403 on this worker host; metadata marks files unrestricted, so acquisition workers should verify from another client/network before bulk downloads. No raw data were downloaded.
- Provenance: deterministic inventory run via `scripts/provenance.py`, seed 0, versions/inputs/outputs in `analysis/runs/`. Next: `data-access`, `data-geno`, `data-meta` should verify downloads, hashes and joins.
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
## 2026-10-08 04:19Z — lean-setup — Lean 4 + Mathlib project bootstrap
- **Task:** lean-setup
- **Agent:** agent-11
- **What was done:** Completed bootstrap of Lean 4 + Mathlib project with trajectory replication framework
- **Achievements:**
  - Created lakefile.lean with Mathlib dependency pinned to master
  - Implemented TrajectorySelection module structure with 4 submodules
  - Defined foundational types (PopulationId, TimePoint, AlleleFrequency, LocusId)
  - Created Trajectory data structure with frequency operations
  - Implemented direction analysis for replicated changes
  - Created frequency path analysis and concordance scoring framework
  - Documented axioms and gaps in lean/AXIOMS.md
- **Validation:** Project builds successfully with `lake build`, all imports resolve correctly
- **Next:** Ready for lean-defs task focusing on core definitions
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
## 2026-10-08 04:07Z — repro-harness — Provenance and run-manifest harness completed
- Implemented comprehensive provenance and run-manifest system for reproducible research
- Created scripts/provenance.py with CLI and library interfaces
- Developed docs/PROVENANCE.md with comprehensive documentation and usage guide
- Added scripts/provenance_example.py for integration examples
- System records: git commit, command line, software versions, input file hashes, parameters, random seed, outputs
- Writes to analysis/runs/<timestamp>-<slug>/manifest.json structure
- Validated system with demo runs and file verification
- Committed detailed implementation and documentation
- Key finding: System provides robust reproducibility tracking for trajectory-replication project
- Next: Integrate provenance recording into all analysis pipeline steps
## 2026-10-08 04:05Z — recon-akbari — Comprehensive Akbari methodology reconstruction completed
- Read Akbari et al. 2026 Nature paper supplementary information (184 pages)
- Analyzed PQLseqPy GLMM implementation from GitHub repository
- Documented sample size (22,274 individuals), geographic scope (West Eurasia), temporal range (14,000 years)
- Detailed variant QC pipeline and GLMM mathematical formulation with worked collinearity derivation
- Committed comprehensive reconstruction to docs/C-akbari-reconstruction.md with seed tracking note
- Found that provenance.py script exists but is in different worktree; noted manual seed tracking approach
- Next: downstream tasks (sim-glmm, lean-collinearity, sim-scaffold) can now build on this foundation
## YYYY-MM-DD HH:MMZ — <task-id> — <one-line>
- what was done
- what was found
- blockers / next
```
## Initial Setup — Repository baseline
- WorksGood task graph initialized
- Research mission defined in `docs/MISSION.md`
- Project structure established with lean, sim, analysis, and docs directories

## 2026-10-08 — sim-glmm — Corrected PQL/LOCO comparator
- Replaced an earlier prototype's individual-row LOCO and incorrectly centered GRM with marker-chromosome LOCO, per-locus binomial PQL/REML variance fitting and Wald Z/p.
- Seeded structured-null, strong single-locus sweep and fixed-variance published PQLseqPy fixture checks are recorded in `sim/test_akbari_glmm.py`; see `docs/F-reproducible-code.md` for commands, assumptions and numerical comparison.
- No locus-by-locus paper summary validation without matched counts, QC and GRM; demographic replacement calibration remains for `sim-sweeps`. Provenance manifest under `analysis/runs/`.

## 2026-10-08 05:02Z — lean-setup — Reproducible Lean bootstrap repair
- Corrected the earlier entry's “pinned to master” claim: `lean/lakefile.lean` and the repository-root Lake entry now require Mathlib `v4.34.1` (locked to `d13f23b...`), with checked-in Lake manifests and Lean 4.34.1 toolchain files.
- `lake build` passed from both `lean/` and the repository root; `lake env lean Test.lean` resolved the imports in both layouts. No project `sorry`, `axiom`, or `native_decide` beyond source comments; `lean/AXIOMS.md` now documents a per-use policy and flags inherited statistical placeholders as non-results.
- Build provenance with seed 0 (no randomness), versions, locked inputs and output hashes was verified under `analysis/runs/20261008-0500-lean-setup-build/` and `analysis/runs/20261008-0502-lean-setup-root-build/`. Next: `lean-defs` can replace placeholder scientific definitions.
