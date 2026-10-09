# Ancient-DNA prioritisation — third-pass results, 2026-10-08

## Headline (incremental, first results)

[observation] **SLC24A5 is NOT recovered high (rank 9,161/10,641, percentile 13.92); TYR narrowly misses the top decile (rank 1,120, percentile 89.48).** Six of eight preregistered control regions reach the top decile. LCT/MCM6 ranks 1, SLC45A2 6, FADS1/FADS2 46, OCA2/HERC2 105, TLR1/TLR6 114 and the AMY1 SNP-proxy region 698. This is a mixed ranking result, not unanimous recovery or a selection discovery.

[observation] The preregistration was persisted at 23:55:38Z before either scorer attempt; successful scoring started at 23:58:54Z. The failed first attempt stopped on an empty pgen-matched path and produced no ranking. The dated implementation clarification and timestamp correction are retained in docs/prioritization-prereg.md. Initial results were written here before drift simulation and bootstrap completion; those results will be appended without replacing this section.

## What was ranked

[observation] data/derived/priority_loci.tsv contains 10,641 supported variants from a 10,642-variant bounded, fixed-250kb-window panel plus predefined control anchors, observed in 2,265 exact-matched ancient people. Every row exposes the directional LLR, magnitude, graded consistency, N_eff-discounted breadth, shape/monotonicity/smoothness and measured composition/carriage components. Variants are not independent causal loci. The genotype cache and sample/panel metadata are in data/derived/prioritization/; fig/prioritization_validation.png is produced after validation.

[inference] This is a PRIORITISER for modern association/phenotype follow-up, not a selection significance study. False positives are expected and acceptable. Neither per-locus selection probabilities nor an FDR are claimed. The frozen null is the previously held-out-calibrated **empirical common-slope variance law**, NOT a generative drift+jump/path likelihood. Shape is a separate descriptive modifier, not a reconstructed selection history. No scored-locus selection amplitude, donor or drift-strength parameter is fitted to optimize the score. Carried-in alleles are down-weighted, never excluded on that basis.

[observation] Source cache is Akbari-ascertained: all QC |X|>=3.16 variants and only 150,000 sampled lower-|X| variants. One closest-midpoint variant per genomic window and predeclared control anchors were selected without our scores, but this does not remove cache ascertainment. Positive-control anchor supplementation favours coverage. Percentiles and Akbari coverage below refer to this panel ONLY, not 9.7 million summary variants or the whole genome.

## Positive controls and precision/recall

[observation] Percentile is 100*(1-(rank-1)/10,641); region rank is the best scored variant within the preregistered GRCh37 interval, not necessarily the causal allele.

| [observation] Control | Best rank | Percentile | Covered variants | Best variant |
|---|---:|---:|---:|---|
| LCT/MCM6 | 1 | 100.00 | 4 | 2_136608646_G_A |
| SLC24A5 | 9161 | 13.92 | 2 | 15_48426484_A_G |
| SLC45A2 | 6 | 99.95 | 3 | 5_33951693_C_G |
| TYR | 1120 | 89.48 | 3 | 11_88874969_T_C |
| OCA2/HERC2 | 105 | 99.02 | 4 | 15_28579961_C_T |
| FADS1/FADS2 | 46 | 99.58 | 3 | 11_61570783_T_C |
| TLR1/TLR6 | 114 | 98.94 | 3 | 4_38625501_G_C |
| AMY1 SNP proxy | 698 | 93.45 | 5 | 1_104012019_C_CA |

[observation] Predefined control variants occupy 27/10,641=0.2537% of the panel. The top 100 contain four control-region variants (4%) spanning three of eight regions (37.5% recall), approximately 15.76-fold above this panel prevalence. Top 10 contain two (20%), spanning two regions. Top decile contains 13 control variants (1.22%), spanning six regions (75% recall).

| [observation] Top N | Control variants | Precision proxy | Regions recalled /8 | Recall |
|---:|---:|---:|---:|---:|
| 10 | 2 | 20.0% | 2 | 25.0% |
| 25 | 2 | 8.0% | 2 | 25.0% |
| 50 | 4 | 8.0% | 3 | 37.5% |
| 100 | 4 | 4.0% | 3 | 37.5% |
| 250 | 8 | 3.2% | 5 | 62.5% |
| 500 | 10 | 2.0% | 5 | 62.5% |
| 1065 (decile) | 13 | 1.22% | 6 | 75.0% |

[inference] **Control precision is a sparse-label proxy, not actual precision for selected loci or modern trait associations.** Unknown selected variants outside these eight regions are not validated negatives; multiple variants in one region are not independent true positives. A claim of 96% actual false positives from the 4% top-100 control precision would be incorrect. True downstream phenotype precision cannot be estimated without downstream labels. Missing SLC24A5 recovery limits sensitivity to known selection forms despite overall enrichment.

## Pending at this incremental write

[observation] Independent random-walk/jump-only tail costs, pure-replacement shape behaviour, 30 bootstrap replicates, carriage-sensitivity ranks and exact Akbari disagreement identities are running or being summarized. No pending number is represented as a completed result.

## Completed validation — appended 2026-10-09T00:02Z

[observation] All three null mechanisms completed 300/300 supported simulated loci; all-panel node-stratified bootstrap completed 30 replicates. Frozen pgen tau²=.0001533793842 and directional amplitude a=.0092278602621 AF/kyr use 389 OLD TRAIN loci, not new candidates. The LLR uses a Gaussian approximation with the calibrated variance; the prior 3.85% gate used a finite-training Student-t threshold, so that gate is NOT an independent validation of this score or shape model.

### Drift-only follow-up cost, not an FDR gate

[observation] Simulation used actual sample dates/path nodes, panel-member frequency baselines and measured between-path frequency differences for jumps. Brownian diffusion .00948676/kyr was normalized to the calibrated COMMON-slope variance using actual-date contrast weights, not fitted to rankings. Jump directions were shared across paths within a replicate, with independently chosen jump times: a deliberately coherent neutral replacement stress case. Binomial dosage sampling is an approximation to imputation and ancient genotyping; these simulations are not calibrated Wright–Fisher histories or a reconstructed Ne.

[inference] Each reported rate means P(a separately simulated neutral locus exceeds the observed Nth-ranked score), i.e. whether inserting ONE synthetic neutral candidate would put it in the observed top N. It is neither a genome-wide FDR nor the fraction of top N that is actually neutral. Tail cost changes with the mechanism and composition of the candidate universe; no arbitrary pass/fail gate is applied.

| [observation] N | Brownian drift only | Drift + replacement jumps | Pure replacement |
|---:|---:|---:|---:|
| 10 | 0/300=0% | 1/300=0.33% | 0/300=0% |
| 25 | 0/300=0% | 6/300=2.00% | 2/300=0.67% |
| 50 | 0/300=0% | 20/300=6.67% | 8/300=2.67% |
| 100 | 0/300=0% | 32/300=10.67% | 18/300=6.00% |
| 250 | 2/300=0.67% | 79/300=26.33% | 62/300=20.67% |
| 500 | 10/300=3.33% | 126/300=42.00% | 103/300=34.33% |
| 1065 | 15/300=5.00% | 165/300=55.00% | 143/300=47.67% |

[observation] At top100 the exact 95% Monte Carlo intervals are 0–1.22% for Brownian drift, 7.41–14.72% for drift+jump and 3.59–9.32% for pure replacement. Zero simulated exceedances is not zero population risk. The corresponding observed score cutoff is .142815; score units have no selection-probability interpretation.

| [observation] Mechanism | Median score | 90th | 95th | 99th | Maximum |
|---|---:|---:|---:|---:|---:|
| Brownian drift | -3.581 | -2.266 | -1.886 | -.802 | -.171 |
| Drift + jumps | -1.730 | .166 | .906 | 2.303 | 4.790 |
| Pure replacement | -1.996 | -.186 | .229 | 1.198 | 1.785 |

[inference] The replacement cost is substantial but acceptable for an exploratory prioritiser if downstream association work can absorb it. The top decile's 55% drift+jump threshold cost argues against treating all 1,065 entries as strongly directional candidates. It does NOT justify excluding carried-in loci or claiming they were not selected elsewhere/earlier.

### Shape: required term present, but biological separation is weak

[observation] The software shape diagnostic assigns 1 to an exact time-proportional ramp and 0 to a noiseless single step; its regression test passes. With sampling noise, pure replacement has shape median .267 (10–90% .048–.377), compared with Brownian drift .258 (.063–.363) and drift+jump .276 (.064–.367). Observed top100 median shape is .297 versus .266 panel-wide. High-ranked non-Akbari-significant chr6:28,875,501 has shape .099 despite rank7; chr1:36,381,622 has shape .036 despite rank10.

[inference] **The shape modifier does NOT convincingly separate noisy replacement from drift in this panel.** A bounded log penalty can be overwhelmed by magnitude evidence; sampling jitter and only three/four nodes limit resolving power. Consequently the ranking is enriched for some known selected regions, but is NOT validated as a strong classifier of sustained ramps versus discrete jumps. This is a partial negative result and the main methodological limit to address before stronger mechanistic language. No score weight was retuned after seeing this limitation.

### Bootstrap stability

[observation] Thirty within-path, within-node individual resamples recomputed ALL 10,642 input-variant scores with the frozen null, amplitude and background relatedness. Mean whole-panel rank Spearman=.8004. Original top-N membership is only moderately stable; per-locus probabilities are in data/derived/prioritization/bootstrap_membership.tsv.

| [observation] N | Mean original entries retained | Mean retained fraction | Mean Jaccard |
|---:|---:|---:|---:|
| 10 | 6.83 | 68.33% | .537 |
| 25 | 12.53 | 50.13% | .340 |
| 50 | 22.30 | 44.60% | .289 |
| 100 | 43.73 | 43.73% | .281 |
| 250 | 112.53 | 45.01% | .291 |
| 500 | 235.97 | 47.19% | .309 |
| 1065 | 560.73 | 52.65% | .357 |

[inference] Sample bootstrap does not assess uncertainty in trajectory inference, dates, imputation or the old TRAIN law; treating rank99 and rank101 as qualitatively different would be unwarranted.

### Carriage weighting sensitivity and breadth

[observation] Against the preregistered carriage floor .25, using floor .5 gives rank Spearman .9958 and top-decile overlap 93.71%; removing the carriage penalty entirely gives .9806 and 86.20%. Both alternative ranks are in priority_loci.tsv. Thus overall order is fairly insensitive to this weighting choice, though about 14% of top-decile membership changes without it. The baseline decomposition remains descriptive, not causal attribution.

[observation] Median supported lineage count is39; the measured-PC kernel participation-ratio N_eff is2.1104, with the same value at the 10th and90th percentiles. Many variants have identical supported sets, so breadth contributes almost a constant to their order.

[inference] This very conservative similarity proxy does not establish 39 independent replicates and offers little discrimination within this panel. It is reported rather than replaced with a more flattering measure after results.

### Agreement AND disagreement with Akbari

[observation] Akbari X is NOT centred; |X|>5.45 is the pi>99% reference threshold, whereas |X|=3.16 is approximately FDR50%, not genome-wide significance. Of108 significant scored variants,73 (67.59%) enter our top decile. The eligible cache contains31,939 significant variants; our coverage of those is only108/31,939=.338%, so 67.59% is explicitly a covered-panel result, NOT genome-wide recall. Whole-panel score versus |X| Spearman is **-.1182**, despite enrichment among the most significant reference variants.

[observation] **992 of1,065 top-decile variants (93.15%) are NOT Akbari-significant.** Exact identities and full components for ALL992 are in data/derived/prioritization/top_decile_not_akbari_significant.tsv. Highest-ranked examples:

| [observation] Variant | Our rank | Akbari absolute X | Our score |
|---|---:|---:|---:|
| 6_28875501_T_C | 7 | 5.190 | 4.753 |
| 1_36381622_C_T | 10 | 4.731 | 4.139 |
| 11_14374308_G_T | 17 | 3.439 | 1.845 |
| 6_25875567_CTGTAT_C | 19 | 4.237 | 1.777 |
| 10_113634916_A_AC | 20 | 1.636 | 1.763 |
| 7_76373970_C_G | 21 | 5.057 | 1.736 |
| 7_100125975_A_G | 22 | 3.791 | 1.610 |
| 18_59876495_A_G | 27 | .671 | 1.397 |
| 11_871030_G_A | 28 | 5.044 | 1.278 |
| 11_47375553_G_A | 30 | 3.464 | 1.147 |

[inference] These are the prioritisation contribution: extra follow-up candidates, NOT demonstrated Akbari false negatives or proven demographic over-correction. Some are close to the Akbari threshold, others have very low reference X. Weak shape at the very highest extra candidates must remain visible rather than hidden by their total rank.

### AMY and SLC24A5 interpretation

[observation] AMY's rank698 is a local SNP/indel proxy, not AMY1 CN selection. The prior CN audit supplies533 coarse-category ancient records but only one adequately supported continuity-constrained CN path; its CN results are preserved in docs/amy-mixture-null-results.md. No SNP tau² or likelihood was applied to CN.

[observation] The included SLC24A5 anchor has ALT AF=.01876, common beta=.000057/kyr and corrected SE=.002664; the observation is nearly fixed toward the opposite allele in this ancient matched panel.

[inference] Historical selection before the observed within-path interval can be missed by a sustained-change prioritiser; that is a scope/sensitivity limitation, not evidence that SLC24A5 selection failed or that demographic carriage disproves selection. TYR's best regional proxy lies just below the decile cutoff, but its causal relevance was not independently assayed here.

### Modern traits and practical recommendation

[observation] A filename inventory of data/raw/ found no modern GWAS catalog, locus–trait association table or phenotype scan. No valid offline modern trait enrichment was computed; gene names and known-selection status were not substituted for phenotype association labels.

[inference] **Recommend a top100 exploratory follow-up pilot**, not a definitive selection panel. Its measured control-region precision proxy is4%, recall3/8 regions, and simulation threshold cost is10.67% under coherent drift+jumps (6% under pure replacement, 0/300 under independent Brownian drift). Only43.7% of original top100 membership persists on average in sample bootstraps, so supply ranks, components and bootstrap membership to the downstream study rather than a hard biological cutoff. If broader coverage is inexpensive, top250 covers5/8 controls but incurs26.3% drift+jump threshold cost.

[inference] **Actual expected modern-association/phenotype precision is UNKNOWN.** Four percent is an observed positive-control proxy, not a promise of4% downstream associations. Estimating the latter requires modern association outcomes; this study cannot honestly supply that number. Downstream testing should filter candidates; ancient ranking alone is evidence for prioritisation, NOT a discovery claim, an FDR estimate, independent selection replication, causal phenotype mechanism or a proof that jump-shaped candidates are unselected.

## Validation and reproducibility

[observation] Commands: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 .venv/bin/python analysis/prioritization_score.py (initial failed-empty-path attempt); the same command with --reuse completes ranking,300 simulations per mechanism and30 bootstraps. A subsequent --reuse rerun adds exact cached-index/pvar-allele verification without altering the registered score. Source-named prioritisation tests pass7/7; measured-mixture suites9/9; encoding/path/process suite11/11. Actual AADR-versus-pgen orientation check passes (ALT correlation+.8103). Provenance uses scripts/provenance.py with git probes explicitly disabled, as required. Run manifests under analysis/runs/ record inputs, seeds, source/output hashes and immutable-large-input sizes.

[observation] No git command was run; no file was staged by this worker. Existing repository staging was not queried. Only analysis/, docs/, fig/ and data/derived/ were written; raw data remained immutable. Existing docs were extended with marked appended sections, not rewritten.

## Final artifact verification — appended 2026-10-09T00:05Z

[observation] Verified rerun completed: all10,642 selected cached pgen indices agree with exact pvar CHROM/POS/REF/ALT IDs; supported ranking is unchanged. An independent full-cache score recomputation reproduces every ranked ID and component; automated artifact checks confirm contiguous sorted ranks, support/shape/N_eff bounds, all simulated exceedance counts,30 bootstrap rows and exactly992 top-decile nonsignificant identities. Evidence: analysis/prioritization-artifact-validation.log; source tests and orientation outputs are in analysis/prioritization-*.log. Scoped source/new-doc whitespace/conflict-marker checks pass without git. Reproducible genotype cache is excluded from staging via its local .gitignore; it is not a required source deliverable.
