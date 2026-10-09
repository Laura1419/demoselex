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

## 2026-10-09 00:03Z — prioritization-third-pass (append-only session entry)

[observation] Corrected objective with an appended MISSION amendment; persisted preregistration BEFORE new rankings. Implemented frozen common-slope directional evidence plus graded consistency, measured-PC N_eff breadth, explicit shape and nonzero carriage weighting. Ranked10,641 variants; validated300 drift-only loci per mechanism and30 whole-panel bootstraps. Source tests7/7, previous mixture9/9 and source/process11/11; real ALT orientation passes. All selected cached variant indices were reverified against pvar on rerun.

[observation] SLC24A5 rank9161 and TYR1120 are not recovered in the decile; six/eight controls are. Top100 proxy precision4%, recall3/8; drift+jump threshold cost10.67%, replacement6%; bootstrap retention43.7%. Noisy shape separation is weak. Of108 covered Akbari-significant variants73 are top-decile;992 top-decile nonsignificant exact IDs were emitted. No offline modern trait enrichment computed. Completed results: docs/prioritization-results.md.

[inference] Suitable for an exploratory top100 downstream pilot, not strong ramp-versus-jump discovery. True modern phenotype precision unknown. No old conclusions were rewritten; no git command, staging or raw-data mutation occurred. Provenance manifests use scripts.provenance with git probes disabled; seeds and source/output hashes are recorded under analysis/runs/.

## Appended 2026-10-09 — block-survey — exploratory candidate survey

[observation] Registered first at 02:18:36Z, downloaded/CRC-verified UCSC hg19 refGene with URL/date/SHA256, built 10,615 ranked tiles from 181,933 scored variants (723 strong). Top-100/500/1000 strong counts 86/188/244. Extracted all 479 supplement trajectory tags and checked coordinates/rounded X against local GRCh37 cache; no independent LD intervals or 347-pruned identities available. Table S5 is polygenic context, not locus mapping.

[observation] Two approved post-registration corrections are explicitly appended: cache evidence reveals the scored moderate band was never sampled, reducing non-strong novelty 14/312/756 to no-cache-evidence 0/13/37; transcript interval unions replace min/max symbol bounds, changing 1,367 annotations / 3 initially displayed novel tiles. Frozen ranks/scores/MHC/scored labels remain unchanged. Tests audit real source gaps, positional calls, tile/threshold rules, cache adjudication and artifact provenance.

[inference] False positives are expected, not hidden: non-strong follow-up counts and saved 32/300=10.67% single-locus drift+jump cost are reported; true linked-block FP cost remains uncomputed. Earlier 43.73% variant membership retention and carriage blindness are carried forward. No discovery or modern trait claim. No git command, staging or commit was performed; parent owns integration. Next: independent modern-phenotype evidence, linked-block null/stability; do not infer novel selection from proximity.

## 2026-10-09 02:59Z — two-channel analysis — appended session record

[observation] Pre-registered both channels, additive total, correlated/nested ordered-series null and inclusive ranking before computations in docs/two-channel-prereg.md. Processed all 538,262 scored variants, with 149,994 controls (six fewer than quoted 150k), 10,634 250kb blocks and 626 control-only blocks. Carriage-class 75,876 (control pass 6.2022%), within-class 60,186 (3.4548%), both-large 29,436 (.3127%). Robust ordering 93 variants/eight any-hit blocks, all both-large/strong (control pass 0%); four ordered maximum-score block representatives (control-block pass 0%). Between evidence is retained, not nuisance-zeroed.

[observation] AMY channels .151399 within (P=.2055), .207252 between, additive .358650 copies/kyr; pooled observed .218218 requires the explicitly reported within-information scaling. Seven-component ordering rho=.7143, weighted C=.5825, correlated drift P=.1072/.1958/.3534, conservative P=.3534: a positive CN/carriage motivating example, but not a significant ordering detection under this null. Independent release stream confirms 2,204 AMY SNPs/max|X|=3.794956889, below 5.45. PEL/Pima/Surui means reproduce; PEL H5=.44706 not .46, Europe-five-pop=.17097 not .18. All details/new conclusions appended in docs/two-channel-results.md.

[inference] Large carriage plus robust ordering is candidate selection evidence, not dismissal; no robust small-within/ordered-carriage class was detected (matching control pass 0%). False positives remain acceptable for phenotype follow-up, with control fractions stated for every screen. Same-genotype Akbari reuse, 43.7% prior retention and unavailable calibrated LD-honest block null remain explicit limits. Five new tests, seven existing tests and full artifact validator pass. Seeded runs and validation have script-native provenance with git disabled; no git command/staging/commit, no raw modification. Parent owns review/commit.

## 2026-10-09 — cn-as-snp — bounded same-cohort copy-number counterfactual

[observation] Preregistered before fitting; supervisor approved the explicit ordering gate and AMY-matched primary cohort rather than absolute full-scan reproduction. Matched225/288 ancients; built individual-level genome-wide LOCO GRMs from1,133 filtered systematic SNPs. Canonical controls beat the released lead AMY SNP by8.75–13.75x in |Z|. Five synthetic haplotype loci and continuous AMY1 were fitted with/without GRM, with199-attempt ancestry-conditioned and unrestricted nulls. Primary rerun yielded identical biological results. Seven new tests, eight existing GLMM regressions and exact artifact/GRM audit passed. No git command, staging or raw edit.

[inference] H1^a/H5 are stronger than the released lead SNP but not above transported X cutoffs; LCT-anchored extrapolations11.493/7.556 clear5.45 only as uncertain extrapolations. No GRM-absorption collapse; conditional P=.165/.385 for these two, unrestricted P=.005/.010. Original-paper verdict unresolved, not an exact reproduction or selection claim. H9 requires an explicitly preserved convergence-failure/longer-iteration diagnostic caveat. Full evidence: docs/cn-as-snp-results.md and data/derived/cn_as_snp/.

## 2026-10-09 — NEW amy-binarized analysis session

[observation] Preregistered primary A=H1^a versus all longer states and <=2/<=3-copy sensitivities before results. Reconstructed 225/288 matched ancients, raw-genotype-verified three genome-wide LOCO GRMs, refitted controls and all cuts on/off, completed 2,388 shuffle fits twice deterministically. Seven source-named tests and eight GLMM regressions pass; outputs and provenance under analysis/ and data/derived/amy_binarized/. No git commands or staging; raw data unchanged.

[observation] Primary B=379/450=.84222; Z on/off=-4.1163/-4.2456, direct |X-equiv|2.7261/2.8117; unrestricted P=.005, ancestry P=.13. Within=.04858, between=.03282, total index=.08140 B-frequency/kyr; within MDE80~.14576/kyr. Released band contains356,335 variants versus31,939 strict hits, with9,285 band-only and723 strong 250kb blocks. Released maximum16.6433 on chr5 exceeds canonical LCT16.5232.

[inference] Pooled binary direction is visible and mostly carriage (94.6% pooled contribution); the sparse GRM does not remove it. Ancestry-centred correction leaves a nonsignificant, underpowered within estimate. Primary binary is exactly the prior H1^a complement, so pooling cannot improve on that test's power. Modern B is high in several agricultural populations but saturates in Mbuti too; diet-specific replication depends on the cut. Full results append in docs/amy-binarized-results.md; false-positive and scale caveats retained.

## NEW SESSION — local-ancestry Stage 1, 2026-10-09T17:03Z

[observation] Wrote preregistration before genotypes; implemented exact-label source panel, relative/duplicate exclusions, held-out split, 22,986-AIM panel, 15-state unphased ancestry-pair HMM, forward–backward/Viterbi, independent-marker PCs, nine sensitivity settings and fail-closed AMY pipeline. Final-code replication reproduces all finite validation statistics; six source-named software tests pass and ALT orientation reproduces r=.8103. Required scientific validation FAILED (IRAN coverage, EHG predominance, proportion stability), so AMY test was stopped. Numeric report: `docs/local-ancestry-results.md`; artifacts: `data/derived/local_ancestry/`; provenance manifests use scripts/provenance.py with ALL git probes disabled. No git commands, staging, commits, package installation or immutable raw modifications performed.

[inference] This is a reproducible negative machinery-validation result, not success at the validated-painter goal. No new AMY null or local detection-floor claim is made. Prior global baseline and raw CN/canonical X facts were rechecked; strict CN+AADR mapping519 differs from supplied522. Further segment-validation/reference-coverage work must be a new preregistered revision.

## APPENDED session — 2026-10-09: full-panel haplotype backgrounds

[observation] Implemented and preregistered all-panel unphased genotype-window clustering with two window sizes, four K values, per-site exclusion audits, structural H1/H2 inventory, sequence-distance joins, detection-floor machinery and fail-closed validation. All16,339 imputed rows have100% usable background dosage coverage. Found a biological exclusion-coordinate defect despite numerical exclusion tests: initial AMY interval missed gene copies. Supervisor approved an explicit post-results annotation-derived repair for AMY and all controls; original outputs remain retained, not interpreted. Corrected known-structure/stability validation remains failed, so no AMY outcome conditioning is interpreted and no canonical recovery is asserted.

[inference] Coverage is solved for this input representation, but stable biologically interpretable local diplotype backgrounds are not demonstrated. No selection, independent-origin, or achieved power claim. New outputs live under analysis/, docs/, data/derived/ only. No git command/staging/commit; parent handles integration.

## APPENDED CONTINUATION — 2026-10-09 19:15Z — local-ancestral-state feasibility lane

[observation] Read the committed method study first, preserved the original preregistration/snapshot, and appended an outcome-blind continuation registration. Added `docs/local-ancestral-state-design.md`: explicit lower-copying-state choice, pinned ELAI 1.01/source revision, missing joint/virtual-site exporter and QC/truth requirements, exact prospective pipeline and honest no-prototype verdict. No new biological conditional statistic was computed.

[observation] Twelve synthetic tests pass, including new fail-closed corrected-AMY exclusion assertions. Fresh chunked released-X audit reproduces 2,204 AMY SNPs, max |X| 3.794956889119673 and zero above 5.45; canonical-locus maxima also match. Headers confirm 52,382,872×16,339; both saved feature manifests have zero forbidden coordinates/indices, and all 22 core transcripts are enclosed.

[observation] Explicit discrepancies: current dated panel metadata has 5,987 unique biological IDs, not ~5,737; saved primary 1 Mb/K=4 silhouette is 0.154051, not blanket <0.10. Prior structure/stability still fail. The ~0.049 slope/~0.146 floor belongs to the primary A≤1-complement long-B outcome, not the A≤2 sensitivity. Neither old global N_eff nor 100% imputed call completeness validates local shared states.

[inference] STOP before inference: new-method structure/stability/coverage/negative-control gates remain unrun/unpassed; local shared-state N_eff and achieved floor unknown. Next requires authentic lower-pair export plus masked-preprocessing/QC and ancestry-matched truth, followed by all preregistered gates—not another k-means or independent-SNP HMM.

[observation] Initial root free space 52G (100% reported used); audit start/end free bytes 58,279,493,632/63,346,241,536. No genotype caches, large matrices, NPY files, PGEN/PVAR copies or scratch created. All writes stay in permitted directories. Provenance uses `scripts/provenance.py` with git probes disabled; no git commands or staging. Parent reviews and commits.

### 2026-10-09 — ELAI lower-state joint posterior patch (GPL-3.0-or-later additions)

Patched the pinned ELAI fork with exact ordered/unordered lower-pair exports and
an export-only missing-emission bridge. Clean build and all checks passed:
15 numerical gates, 4 movement gates, 5 rejection cases, 18 baseline byte
comparisons. Documented upstream --pk1 copy-count scaling and admixture-output
updates, bridge homogeneity/sensitivity, and absence of real aDNA validation in
`docs/elai-patch.md` / fork `README-PATCH.md`. No files staged or commits made;
parent owns independent acceptance and publication. Provenance recorded under
`analysis/runs/` with slug `elai-lower-state-posterior-final`.

## APPENDED peer-review evaluation — 2026-10-09T20:02:50+00:00

[observation] Read all61 physical pages of Akbari et al.2026 Peer Review File using pypdf; verified immutable PDF SHA256. New docs/akbari-peer-review-evaluation.md inventories V0/V2/V3,99 substantive concern rows,264 page/version/speaker-checked verbatim excerpts, five-topic deep dive and all six requested project cross-checks. Independent re-extraction and source/section/table/absence checks pass. No git command, staging or commit performed; index not inspected. Parent reviews and commits.

[inference] Review repeatedly challenged threshold calibration/false positives and explicitly raised time–GRM collinearity, but did not explicitly ask about removal of true ancestry-aligned selection. AMY/copy number/Bolognini/Sudmant/Garrison and the3.16-to5.45 prioritization band are absent from accessible prose. Source-only selection is an explicit negative control for the authors’ transect estimand. Preserve our negative replication result; do not treat the review as adjudicating AMY or proving an over-correction class. Cross-check flags that X-equivalent3.64 is A<=2 sensitivity, whereas97% GRM retention belongs to primary A<=1 (X-equivalent2.726).
