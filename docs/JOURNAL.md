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

## 2026-10-08 05:14Z — lean-collinearity — Time/GRM mode geometry
- Proved `||r||² = ||t||²(1−cos²θ)`, a near-alignment bound, and the hard-projection VIF `1/(1−cos²θ)` for nonzero residual in `lean/TrajectorySelection/Collinearity.lean`.
- The finite one-mode working model retains `α²/(σ²+λ)` information even at perfect alignment. Corrected Section C's earlier claim that finite GLMM fitting exactly projects away `col(K)`; Section D states the linearized covariance assumptions and limits.
- Lean build and focused tests passed without new `sorry`/axioms; `#print axioms` reports only standard Lean foundations. Provenance (seed 0, Lean 4.34.1 and pinned Mathlib) is in `analysis/runs/20261008-0514-lean-collinearity-proof/`. This is a conditional linear-algebra result, not empirical evidence of over-correction.

## 2026-10-08 — sim-glmm — Corrected PQL/LOCO comparator
- Replaced an earlier prototype's individual-row LOCO and incorrectly centered GRM with marker-chromosome LOCO, per-locus binomial PQL/REML variance fitting and Wald Z/p.
- Seeded structured-null, strong single-locus sweep and fixed-variance published PQLseqPy fixture checks are recorded in `sim/test_akbari_glmm.py`; see `docs/F-reproducible-code.md` for commands, assumptions and numerical comparison.
- No locus-by-locus paper summary validation without matched counts, QC and GRM; demographic replacement calibration remains for `sim-sweeps`. Provenance manifest under `analysis/runs/`.

## 2026-10-08 05:13Z — lean-replacement — Sign-null and replication separation
- Added fully proved two-trajectory finite sign model in `lean/TrajectorySelection/Replacement.lean` and regression examples in `lean/TestReplacement.lean`.
- Constant within-source trajectories plus an **additional iid fair-sign noise assumption** give the exchangeable-sign null. Shared positive true drift plus an **additional iid sign-detection bias `q>1/2`** gives `P(C=1)=q²+(1−q)²>1/2`; drift alone cannot justify the sign law.
- `cd lean && lake build` and focused Lean tests pass. `#print axioms` reports only Lean's standard `propext`, `Classical.choice`, `Quot.sound`; no project `sorry` or new axiom. This is a two-trajectory model, not a calibrated demographic scan. Section D records scope and limitations. Seed-0, version-pinned provenance: `analysis/runs/20261008-0514-lean-replacement-compile/manifest.json`.

## 2026-10-08 05:02Z — lean-setup — Reproducible Lean bootstrap repair
- Corrected the earlier entry's “pinned to master” claim: `lean/lakefile.lean` and the repository-root Lake entry now require Mathlib `v4.34.1` (locked to `d13f23b...`), with checked-in Lake manifests and Lean 4.34.1 toolchain files.
- `lake build` passed from both `lean/` and the repository root; `lake env lean Test.lean` resolved the imports in both layouts. No project `sorry`, `axiom`, or `native_decide` beyond source comments; `lean/AXIOMS.md` now documents a per-use policy and flags inherited statistical placeholders as non-results.
- Build provenance with seed 0 (no randomness), versions, locked inputs and output hashes was verified under `analysis/runs/20261008-0500-lean-setup-build/` and `analysis/runs/20261008-0502-lean-setup-root-build/`. Next: `lean-defs` can replace placeholder scientific definitions.

## 2026-10-08 05:06Z — lean-defs — Generic trajectory and concordance definitions
- Added `lean/TrajectorySelection/Defs.lean`: bounded indexed frequency paths, weighted OLS slope with undefined-case handling, exact direction and nonzero-majority concordance, constant-source turnover mixture, and PSD covariance kernel. No new Lean axioms or sorry.
- Added focused Lean examples in `lean/TestDefs.lean`; updated Section D with equations and explicit limits, retaining inherited placeholder warnings.
- Both root and `lean/` builds passed; the only warnings came from inherited scaffolds. Deterministic build provenance uses seed 0 and pinned Lean/Mathlib versions; no randomness. Downstream theorem tasks must state null/independence assumptions.

## 2026-10-08 05:16Z — lean-statistic — Family-discounted majority statistic
- Proved count concordance range and fixed-`n` monotonicity; defined correlation-mass `N_eff` on *distinct inferred families*, proved `[1,m]` range and antitonicity as relatedness rises. Adding a near-copy assigned an existing family leaves `N_eff` exactly unchanged.
- Fair independent-sign null has a folded-binomial mean and variance; Lean proves exact finite-sum formulas and the fold, while Section D documents the central-binomial closed forms and cautions that exchangeability alone is insufficient. A naive raw-label design effect can increase on duplicate; family assignment is an explicit assumption, not inferred by the proof.
- Focused `lean/TestStatistic.lean` checks ties, zero denominator, `n=2` moments and duplication; build provenance uses seed 0, pinned Lean/Mathlib and version/hash manifest in `analysis/runs/`.

## 2026-10-08 05:23Z — sim-stat — Weighted trajectories and correlated meta-analysis
- Implemented `sim/trajectory_stat.py`: pooled Model 0 and within-trajectory count-weighted slopes, binomial SE, raw concordance, preassigned-family collapse, relatedness-discounted GLS, effective family count and heterogeneity. Section D states assumptions and Section F the reproducible entry point.
- Designed seed-31415 positive and negative signals recover the expected signs; a turnover-only fixture produces a positive pooled slope with zero within-trajectory slopes. All 14 targeted simulation/comparator tests pass; `python3 scripts/provenance.py verify analysis/runs/20261008-0523-sim-stat-validation/manifest.json` passes with hashed inputs/outputs and recorded Python/NumPy/pandas/SciPy versions.
- Default distinct-family independence is only a working assumption. Binomial and GLS p-values require matched neutral calibration before any selection claim; downstream `sim-sweeps` and `sim-ctrl-null` own that test. No new Lean axioms or unseeded random draws.

## 2026-10-08 05:28Z — sim-ctrl-null — Matched-control empirical rank
- Added deterministic same-chromosome six-feature caliper/nearest matching, 1 Mb linked-window exclusion and exact observed trajectory/time support; all candidate and control counts use the same fitted statistic. Candidate loci cannot be controls; insufficient or undefined matches fail closed. Conservative +1 empirical P and an auditable midrank percentile are reported.
- Seed-251008 neutral four-stratum holdout gave 128 candidate ranks against 40 matched controls each: mean P 0.4825, uniform KS D 0.1027, fraction P≤.05 0.0625. Focused matcher/statistic tests (10) and provenance verification pass. See Section G and analysis/runs/20261008-0528-sim-ctrl-null-validation/; these toy independent-locus neutral results do not establish calibration under realistic replacement, ascertainment or ancient-DNA QC.
- No Lean changes, new axioms/sorry or unseeded randomness.

## 2026-10-08 05:39Z — lean-collinearity — Finite-mode near-alignment bound (isolated retry)
- Existing hard-projection and finite-mode results had already landed in main. Added the conditional finite bound `I_time/I_isotropic ≤ ε+δ` for residual fraction at most ε and mode-retention fraction at most δ. This distinguishes large GRM variance from geometric alignment; it does not claim finite-variance non-identifiability.
- Focused `(3,1)` example with mode variance 99, noise variance 1 gives at most 11% baseline information. Full library and focused Lean test pass; theorem's axiom print shows only standard Lean foundations. Provenance uses seed 0 (no random draws) and pinned Lean/Mathlib versions.

## 2026-10-08 05:00Z — sim-scaffold — Repaired and validated dated mosaic simulator
- Replaced the earlier final-only, non-donor migration prototype with sampled diploid mosaic Wright–Fisher trajectories, actual donor migration/turnover, and signed selection.
- All six fixed-seed configs now emit aligned dated genotype/sample matrices, per-trajectory locus truth, and `scripts/provenance.py` manifests; end-to-end tests pass. SLiM was unavailable.
- Pure replacement has zero selection truth despite a locus-frequency rise; opposing selection has negative truth in pop2. These are scenario controls, not power estimates.
- Next: sim-sweeps should vary seeds, migration and turnover to measure calibration, power and model misspecification.

## 2026-10-08 05:59Z — sim-scaffold — Retained scaffold revalidated
- Three simulator tests pass, including two end-to-end runs of all six fixed-seed scenarios and checks of dated genotypes, truth, and provenance; the required diff whitespace check passes.
- The retained branch is clean before this journal entry but does not contain advancing `main`; WG completion previously rejected it on ancestry, and project policy reserves branch integration for the coordinator. No worker merge/rebase/pull was attempted.
- Next: obtain a sanctioned reconciliation of the retained candidate, then retry WG completion without relaxing validation.

## 2026-10-08 20:10Z — independent directional-AMY reanalysis (bounded session)

[observation] Preregistered anchored paths and a held-out Gaussian drift-plus-jump check before new outcomes, with timestamped numerical/feasibility amendments. Corrected load_aadr.py's verified TGENO bit-order error; ten source-named regressions pass, including actual HO/1240K call concordance and full random-walk covariance.

[observation] Rebuild: 42 candidate paths, 35–377 people, 3–5 nodes, 2,151–4,850 years. Ancestry effective diagnostic 15.90, not a genealogical replicate count. Six-path held-out SNP pass fraction 13/208=.0625 (95% interval .0337–.1045), preliminary calibration gate passed; Ne fits at infinity boundary and is not resolved. AMY: 522 numeric AADR matches, 87 carriers in 13 paths, only one long supported CN series, zero qualifying >=3-deme successions. CN process probability remains unidentifiable/un-calibrated.

[inference] The prior capstone overstates falsification; locus-rank scans/finite P floors and bad units did not test this CN argument. Conversely, shared neutral donors can manufacture consistent directions, so fair iid signs require evidence. Evidence cannot currently establish or refute AMY parallel selection. See the three new docs and analysis/deme_rebuild/ outputs.

[observation] No git command/staging/commit, raw-data modification, or Lean edit. RunRecorder git probing explicitly disabled; manifests record code/prereg/metadata/output hashes and raw genotype path/size. Prereg snapshots reconstructed from the append-only text match their recorded pre-run hashes. Two initial implementation failures (integer output casting, boundary zero variance) were repaired before successful fitting; no rank transform forced calibration.

[inference] Next: parent review, rerun affected old loader-dependent results, establish donor-resolved continuity and CN haplotype/mutation calibration with more independent long series (including the claimed Andes replication). Do not promote this feasibility result into a selection discovery.

[observation] Before finishing, added separately timestamped pre-run secondary weakest-replicate statistic (Amendment 4), because a pooled sum can be outlier-driven. Held-out directional controls pass **9/208=4.33%** (95% interval **2.00–8.05%**); preliminary gate passes without a post-hoc control-rank transform. Conditional power is **7.53/15.27/42.77/90.22%** at .25/.5/1/2 per-deme SD. Eleven regression tests now pass, including the one-outlier contrast. This secondary result is explicitly exploratory relative to the original preregistration; no biological locus interpretation is enabled.

## 2026-10-08 — second-pass measured-mixture null audit

[observation] New pre-run procedure/results: docs/mixture-null-prereg.md, docs/amy-mixture-null-results.md, docs/genomewide-mixture-null.md. Implemented measured-label CN/SNP residual tools and source-named regressions; recorded deterministic conditional power, pure-composition specificity, per-locus held-out controls, and loss of power as nuisance partitions proliferate. No git command, staging or raw-data edit by this worker.

[observation] AMY has533 dated labeled individuals, EIGHT labels (seven listed categories plus singleton CHG); residual beta .151 copies/kyr, P=.206. Conditional null4.1%, .45 copies/kyr power85%; singleton correction unidentifiable. Strict CN replication infeasible:87 carriers in13 rebuilt paths, only ONE qualifying supported series. Genome-wide raw strict-deme controls fail8.10%/7.95%; one independently TRAIN-fitted variance parameter yields held-out6.08%/3.85% (intervals include5%), median floors .0309/.0260 frequency/kyr. New regression suites5+4 pass;11 existing tests and real AADR/pgen orientation check pass.

[inference] Compatible with coarse composition, not proof AMY was purely demographic or that s≈.02 is false. Corrections should freeze measured source/ancestry predictions, test within-lineage excess against independent measured variation, and publish uncertainty/power. Scope is a bounded genome-wide-distributed control scan, not full discovery or fractional-ancestry/CN evolutionary modeling. Future validation needs untouched LD-block controls and independently measured longitudinal source ancestry/CN.
