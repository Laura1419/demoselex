# Local ancestral sharing-state frequency change: prospective registration

## NEW REGISTRATION — 2026-10-09 (written after method study, before new real-data statistics)

[observation] This is a **conditional design, not permission to interpret an unvalidated AMY fit**. `docs/local-ancestral-state-method.md` supplies the method choice and constraints. No real-data ancestral-state inference or conditional outcome statistic has been run for this registration. The registration will be snapshotted and hashed by a preflight audit before inspecting any new biological statistic. Changes after this point must be APPENDED, dated, and distinguish outcome-blind amendments from post-result tuning.

[inference] Chosen design: authentic **ELAI/fastPHASE lower-haplotype-state inference** with unknown phase marginalized, followed by an **ASAFE-type diploid convolution likelihood with time-varying state-specific frequencies**. Lower states operationalize local haplotype sharing; they are NOT genealogically dated MRCAs. ELAI's standard `--ps2` upper-pair output is not the required lower-pair output. **Implementation is STOPPED pending that interface, uncertainty/QC information and calibrated locus-gap bridging.** No upper-ancestry regression, SNP k-means or hand-built window HMM may be substituted silently.

## 1. Inputs, cohort, populations and outcomes

[speculation] Primary dictionary is outcome-blind, common to ALL 16,339 imputed panel individuals, with no date or AMY structure used for fitting. Use the single-upper-layer linked-model limit (S=1) to define sharing rather than asserting five named sources. Baseline lower K=5 follows the manual's 5×S heuristic; sensitivities K=3 and 8, seeds 20261009/20261010/20261011. This is a heuristic to be tested, not evidence that five ancestral states exist. If S=1 or lower posterior export is unsupported by the selected executable, STOP; do not fabricate support from the mathematical paper. Run 50 EM steps and check convergence by an additional 20 steps (relative log-likelihood change <10^-4 per observed allele); permit only the preregistered grid, not post-AMY K selection.

[speculation] Deduplicate biological individuals and exclude first-/second-degree relatives from dictionary fitting and frequency inference, retaining explicitly flagged assignments for excluded samples. Fit participant-balanced weights (one biological individual one unit), with ancient/modern and technology strata audited, not 503 modern sequences treated as ancient references. Retain all eligible ancient geographic variation, but modern individuals contribute no ancient temporal observation. Preserve source/data IDs and uncertainty, and report the all-panel and post-QC denominators separately.

[speculation] Primary outcome is measured AMY short-vs-long **B=1 if A>1** (primary established definition); sensitivity B=1 if A>2. Count two unordered structural outcomes per diploid individual; no H1-to-SNP-chromosome assignment is asserted. The ~225 dated measured structural individuals are the maximum available overlapping outcome cohort, not 5,737 new structural measurements. Imputed SNP time-series can separately use all uniquely dated eligible ancient carriers (~5,737 indicated by the operator), but must not stand in for measured AMY structure. No CN-as-SNP imputation without an independently validated structural predictor.

[speculation] For cross-population comparison, use **AADR recorded country** as the fixed, mutually exclusive primary population partition (missing country excluded from replication, not filled from outcomes). This is an auditable geographic yardstick, not a claim that a country is a genetic population. Freeze a crosswalk of IDs, dates, countries, technology, publications and kinship before outcomes; report temporal turnover within countries as residual confounding. Published AADR ancestry/group labels are validation information, NOT the fitted local sharing states. Geographic and ancestry-component partitions must not be swapped after seeing AMY results.

## 2. Locus exclusion and state assignment

[observation] Authoritative exclusion E is the inclusive **GRCh37 chr1:103,915,000–104,327,000** interval. It is the corrected core-plus-buffer exclusion, not an additional ±25 kb around those bounds. AMY1 transcripts end through 104,301,314; do not reuse the former 104,255,000 bound. All focal/core-overlapping gene transcripts must be enclosed using the existing hg19 refGene annotation.

[speculation] Baseline features are quality-checked biallelic SNPs in [103415000,103914999] and [104327001,104827000]: 500 kb each flank. Sensitivities use 250 kb and 1 Mb each flank. Model linked markers, do not randomly project/average all dosages or greedily remove LD as if it were nuisance. MAF≥.01 and INFO≥.8 are **proposed QC minima**, motivated by the ancient benchmark, with genotype posterior≥.90 for hard-call-only input; outcome-blind ancestry/date/technology-stratified error benchmarking remains mandatory. Missing genotypes stay missing. If only mean dosage is available, do not reconstruct GP by rounding. A calibrated hard-call lane may proceed only with downsampled ancestry-matched truth and error propagation that passes gates below.

[speculation] Define the state at **x*=104121000**, the midpoint of E, as latent lower-cluster membership of each chromosome at a **missing-emission virtual site**. No tested genotype, copy number, structural distance/label, date or population is an emission/training input. Estimate lower-state unordered pair posterior

`q_i(a,b) = P({Z_i1(x*),Z_i2(x*)}={a,b} | flank genotypes, fitted linked model)`, `a<=b`.

[speculation] Transitions across the 412 kb gap must come from the real linked model and be calibrated against genetic-map distance/known mosaics, not a new ad hoc physical-window switch penalty. Upper-layer averages or independent products of chromosome marginals cannot replace q. A lower cluster's local allele-probability vector/dictionary and same-run transition parameters define its identity. Across populations use exactly the SAME dictionary. Across seeds/windows/K, align states by held-out flank emission similarity (Hungarian alignment at fixed K); for K changes compare co-sharing matrices rather than treating split labels as identical ancestors.

[speculation] A diplotype is usable only if the maximum unordered-pair posterior is ≥.90 AND simulated central-gap assignment accuracy for its support stratum passes ≥.90. Hard pair calls are for confusion/coverage diagnostics; frequencies marginalize q. Neither chromosome is claimed paternal/maternal or attached to H1/H2. Report left-only versus right-only bridge disagreement and central posterior entropy; when a gap is not identifiable, abstain. Very fine states that fix the outcome invariant cannot identify a within-state change.

[speculation] **Required runtime assertions:** all feature positions are strictly outside E, all core transcripts enclosed, tested structural/SNP outcomes absent from every state-fitting input, and the union of their variant indices intersects the feature set in exactly **0** elements. Holding outside-E input and seed fixed, replacing/permuting ALL E genotypes and AMY outcomes must give **bit-identical q**. This is computational exclusion only. If phase/imputation was originally informed by E, retrieve/recompute masked preprocessing or demonstrate error-calibrated independence with truth; otherwise flag the leakage gate NOT PASSED. No claim of statistical independence of linked flank alleles follows from coordinate exclusion.

## 3. Statistic: change IN frequency, not ancestry fraction

[inference] Let `t=-BP/1000` (forward-time kyr), population r, state a. With `p_ar(t)=logistic(alpha_ar+beta_ar*(t-t0))`, pair-specific genotype probabilities are the coefficients of

`[(1-p_ar(t))+p_ar(t) z] [(1-p_br(t))+p_br(t) z]`.

[inference] For unphased genotype/structural B count `g_i ∈ {0,1,2}`,

`ell = sum_i log(sum_{a<=b} q_i(a,b) P(g_i|a,b,r_i,t_i))`.

[inference] For genuinely available genotype read likelihoods, replace P(g_i|...) by `sum_g GL_i(g) P(g|...)`; posterior GP must be prior-adjusted or integrated as a calibrated measurement model, **not mislabeled GL**. Integrate age distributions rather than treating heterogeneous calibrated-date uncertainty as zero. The convolution assumes exchangeable alleles conditional on state/population/time and random pairing conditional on the modeled diplotype. Test its misfit against family/inbreeding and quality strata. This is our temporal/soft-assignment extension of ASAFE, not published ASAFE software behavior.

[speculation] Primary estimand is the 2 kyr change

`Delta_ar = p_ar(-4) - p_ar(-6)` (6,000 → 4,000 BP).

[speculation] Report `Delta_ar/2` in frequency/kyr along with endpoint frequencies and a family-/individual-block bootstrap CI. This is distinct from logistic beta's log-odds/kyr. Only interpolate: state/population must have support in BOTH [6,500,5,500] BP and [4,500,3,500] BP and the fit uses dated eligible observations in [6,500,3,500] BP. No moving the endpoints to maximize an effect. Nonlinear/bin-frequency checks are descriptive sensitivities, not replacements for the primary contrast. If AMY sample support is insufficient, the answer is underpowered/no eligible state, not a more convenient interval.

[inference] The genotype convolution has phase-free identifiable and non-identifiable regimes. Require full-rank observed information for the named state contrasts, condition number≤100 after parameter scaling, bounded two-sided profile intervals, and no equally supported allele/state-swapped mode. Full-rank mean ancestry design alone is NOT sufficient. Perform multistart optimization and pure-state/mixed-pair simulations. No arbitrary ridge or prior may create a state-specific sign in a non-identifiable design.

## 4. “Shared”, directional replication and independence

[speculation] Freeze state eligibility from flanks, dates, populations and quality **without AMY g**. For each state/population/endpoint bin compute expected chromosome count `C=sum_i E_q[n_ia]` and allocation Kish support `H=(sum_i E_q[n_ia])^2/sum_i E_q[n_ia]^2` over independent biological individuals. H is a support diagnostic, not an error-corrected effective sample size. Require **C≥20 AND H≥10 independent individuals in EACH endpoint bin**, posterior/bridge gates, and at least 10 informative homo-state diplotypes over the analysis interval or demonstrably identified likelihood by the tests above. Phenotype/purity outcomes are not used for this eligibility decision.

[speculation] A state is **shared** if at least **two distinct countries** meet the same criteria using the common dictionary. No separate population fits followed by label-number matching. Report the eligible state×country matrix, actual per-bin denominators and number m of shared states. Countries contribute contrasts within a state, NOT m independent ancestral replicates. Genetic links between countries sharing the state must be modeled.

[inference] Per state report all signed `Delta_ar`. Define conservative common-direction magnitude

`T_a = max( min_{r in R_a} Delta_ar, min_{r in R_a}(-Delta_ar) )`.

[inference] `T_a>0` exactly when every eligible population's point contrast has the same nonzero sign; its magnitude is the smallest absolute concordant contrast. Inference uses simultaneous bootstrap/null-calibrated intervals for the contrasts, not multiplying marginal sign P values. Report direction discordance as a result, not delete a disagreeing country. Different selection regimes could legitimately cause discordance; absence of concordance is not universal absence of selection.

[speculation] Fit/report contrasts per state and population, with assignment/imputation/age uncertainty propagated by rerunning the dictionary and likelihood in each outer bootstrap/calibration draw. Estimate the covariance of state scores `T_a` under **joint drift-only null simulations** retaining migration, shared lineages, dates, state occupancies, relatives and the imputation/assignment process. Reuse common simulation replicates across states so the off-diagonal covariance is estimable. With `R` the normalized positive-semidefinite null covariance, report

`N_eff = (tr R)^2 / tr(R^2) = (sum lambda)^2 / sum lambda^2`,

and rank/eigenvalues/Monte Carlo uncertainty. This spectral participation ratio counts independent score dimensions, NOT literal known independent ancestors. For aggregate testing use the joint simulated distribution (e.g. max T or number of concordant states), not binomial(m, .5), and never round N_eff into a count of independent P values. Large negative/off-diagonal dependence is retained, not clipped to inflate replication. Zero-variance states are untestable; covariance not estimable → N_eff unavailable and no replication claim. If N_eff<2, no evidence claim of multiple independent shared-state replications even if many countries agree.

[observation] **Current registered run has 0 validated/evaluable shared states; the biological m and N_eff are UNKNOWN, not estimated as zero.** No lower states have been assigned here. The previous seven-global-component N_eff≈2 is a different population-score calculation, not an estimate for these local states. The absence of a numerical local N_eff is itself a STOP condition; quoting seven countries/components would be misleading.

## 5. Null, gates and locked decision rule

[inference] Descriptive null: `Delta_ar=0` for each state/population (composition can change arbitrarily while conditional frequencies are constant). Selection-inference null: **neutral within-state drift plus migration/relatedness and observation error**, with no directional selection; drift permits realized nonzero changes. Time-label permutation alone is NOT calibrated for ancient turnover or date-correlated genotype error. Use separately calibrated demographic simulations with neutral matched loci and independently specified histories. If neutral migration of finer backgrounds within a state produces concordance, it remains part of the null, not evidence for selection.

[speculation] Gates are design thresholds, NOT literature facts, and all are required BEFORE reading AMY contrasts:

| Gate | Required measurement and threshold | Stop condition |
|---|---|---|
| Implementable state/likelihood interface | Authentic lower K×K joint posterior; normalized pair sum=1, correct diploid symmetry; tiny-case enumeration and missing-gap simulations agree; real source/export verified | Only upper dosages/pairs or invented phase available |
| Known structure (a) | Fit dictionary without labels. Use a fixed stratified 70/30 reference calibration/heldout split among the precedent's five broad AADR source brackets, requiring ≥2 heldouts/source. Train label predictor only on the calibration split's state dosages; publish five-way heldout confusion matrix INCLUDING abstentions, macro recall≥.80, each source recall≥.60; hard-diplotype/group NMI≥.10 on heldout mapped references | Missing source, inadequate holdouts or failed metrics. This diagnoses structure, not segment truth |
| Local truth / errors | Held-out ancestry-matched phased mosaics with real sample-like missingness, damage/imputation and central gap: ≥.90 balanced local pair accuracy; posterior .9–1 bin ≥.90 empirical accuracy, absolute calibration discrepancy≤.05; every source/age/technology stratum error≤.10 | Global label agreement without local truth cannot pass |
| Stability (b) | Fixed K seeds/window/extra-EM: aligned hard-diplotype ARI≥.80; median total-variation L1 of aligned q≤.20 and 90th percentile≤.40. K=3/8 vs primary: chromosome-sharing matrices Pearson r≥.80; state merges/splits mapped explicitly. Left/right central predictions pair agreement≥.80 among supported individuals | Do not equate ARI with prior global L1. ANY prespecified sensitivity failure stops; no rescue by selecting a convenient K/window |
| Coverage (c) | Report all 16,339 input, post-QC fitting, posterior-supported assignment, metadata-mapped, dated SNP and measured structural n; usable q for ≥.95 of input and ≥.90 of each heldout/date/technology stratum | Imputed call completeness alone is not usable-state coverage |
| Locus exclusion/negative control (d) | 0 coordinate/index intersections, all AMY core transcripts covered, forbidden-data perturbation bit-identical q, masked-preprocessing/truth leakage check; 999 independently permuted balanced pseudo-variants have empirical false-positive rate≤.065 at .05 | Missing quality/leakage evidence is NOT a pass; controls must not use tested AMY outcomes |
| Identifiability/null | Information/profile tests above; ≥1,000 null replicates per matched design, empirical type-I≤.065 at .05; simulate temporal state-assignment errors and country turnover; report covariance uncertainty | A constant-frequency “mixture-only” null alone is insufficient to claim selection |
| Sharing/power | ≥2 eligible shared states and N_eff≥2 for independent replication claim; each tested Delta has simulated ≥80% power at its stated detection floor with interval coverage≥.93 | No eligible states is a valid negative feasibility result; no outcome analysis after a failed assignment gate |

[speculation] **Decision rule:** any missing/failed implementation, state-validation or exclusion gate ⇒ STOP and no biological conditional frequencies/slopes/P values. Passing assignment but failing shared support ⇒ report no evaluable replication; frequencies only if independently identified and power-labelled. Passing every gate ⇒ report all predeclared state/population contrasts and covariance-calibrated concordance with multiplicity accounted by the joint null; label as prioritization evidence, not proof of selection or recurrent origin. Underpowered nonsignificance ⇒ unresolved, not “composition fully explains AMY.”

## 6. Detection floor, specified before outcomes

[inference] For each Delta, simulate on the frozen dates, q, kinship, missingness, population sizes and calibrated genotype/state errors. Test a two-sided .05 contrast; baseline within-state frequencies {.1,.5,.9}, effect grid increasing by .01 with endpoints restricted to [0,1]; 1,000 simulations/grid value, report power and binomial Monte Carlo CI. The smallest |Delta| with power≥.80 **and lower 95% Monte Carlo bound≥.80** is the achieved contrast-specific floor. Estimate the all-population common-direction T floor jointly, NOT assuming independent country tests. If .80 power cannot be reached for feasible frequency endpoints, state no detectable effect in that design.

[inference] Independent-binomial approximation, for orientation only:

`MDE80_Delta ≈ (1.959964+.841621) sqrt(p(1-p)/(2 n_early)+p(1-p)/(2 n_late))`.

[inference] Even 225 independent diploid outcomes, perfectly balanced 112/113 between endpoints with a single perfectly known state and p=.5, give an **optimistic ~.132 frequency change / .066 per kyr** floor. With 5,737 SNP-outcome individuals balanced 2,868/2,869 the corresponding floor is ~.0262 / .0131 per kyr. These are synthetic design illustrations, not achieved local AMY power, and cannot be reduced by pretending 225 outcomes exist in each state/population. After splitting by state AND country, relatedness and missing dates, the floor rises; genotype uncertainty further worsens it. The previous .145764/kyr global-history HC3 floor is a different model/interval calculation, not a transferable local floor.

[observation] **Current local floor: unavailable (state assignments and null calibration not implemented/validated).** Prototype status: STOPPED BEFORE INFERENCE, not “gates passed because we did not test them.” This registration defines what would be computed only after the missing real-tool inputs and interfaces exist.
