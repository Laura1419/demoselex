# Local ancestry Stage 1 results

## NEW INCREMENT 1 — 2026-10-09T16:59Z

[observation] `docs/local-ancestry-prereg.md` persisted before genotype analysis. Implemented `analysis/local_ancestry.py`: 15 unordered ancestry-pair states, tempered SNP likelihoods, chromosome-reset forward–backward and Viterbi, missing/low-support exclusion, five exact-label sources and held-out individuals. AIM panel and all nine window/penalty settings are saved in `data/derived/local_ancestry/`.

[observation] First run: 22,986 AIMs from 71,552 source-divergent candidate SNPs; 278 structural-haplotype targets mapped to AADR; 19/27 held-out references evaluable. Overall evaluable dominant-label classification is19/19, but both IRAN holdouts lack coverage, EHG median correct-source fraction=.514<.60, and sensitivity fails (1Mb/penalty2 median global L1=.780>.25). Negative control WHG fraction=.958; independent-PC mean R²=.535. Validation therefore FAILED. Final-code replication and carried-forward fact checks pending; no AMY regression/frequency result has been produced.

[inference] **AMY analysis STOPPED under the preregistered rule.** Accurate dominant labels in the covered subset are not sufficient validation of local ancestry proportions/segments. This is a machinery implementation with a negative validation result, not the requested validated painter. It cannot answer whether AMY changes within a fixed ancestry background or revise the prior carriage-dominated interpretation.

## NEW COMPLETED VALIDATION / STOP VERDICT — 2026-10-09T17:02:13.601732+00:00

### Reference panel and target cohort

[observation] Exact-label raw counts, exclusions and eligible BP ranges (years before1950) follow. Assessment=`Pass*`, all annotated family relations excluded, biological-ID duplicates collapsed using 1240K-hit priority. Eligible totals:105 references,78 training/27 held out. Labels are empirical brackets of variation, not assertions of discrete replacements.

| source | label | published_n | eligible_n | published_min_BP | published_max_BP | train_n | holdout_n | min_BP | max_BP |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| WHG | France_Mesolithic | 18 | 13 | 7187 | 13027 | 9 | 4 | 7187 | 13027 |
| WHG | Luxembourg_Loschbour_Mesolithic | 2 | 1 | 8025 | 8025 | 1 | 0 | 8025 | 8025 |
| EHG | Russia_Karelia_Mesolithic | 16 | 16 | 7200 | 8296 | 12 | 4 | 7200 | 8296 |
| FARM | Turkey_N | 68 | 40 | 5615 | 8661 | 30 | 10 | 5615 | 8661 |
| IRAN | Iran_GanjDareh_N | 17 | 7 | 9800 | 10162 | 5 | 2 | 9800 | 10162 |
| STEPPE | Russia_Samara_EBA_Yamnaya | 46 | 28 | 4518 | 5165 | 21 | 7 | 4518 | 5165 |

[observation] Structural archive has288 dated labelled ancient diplotypes;278 unique biological individuals match AADR,278 remain after training-reference exclusion,278 meet genome-wide coverage criteria, and239 meet the primary AMY-vicinity support rule. Exactly383 reference+target individuals were computationally painted; only19/27 holdouts were evaluable. These are NOT the533 copy-number ancients or the previous225-person imputed cohort. The239 supported AMY targets are diagnostic coverage counts, NOT a tested cohort. `target_mapping.tsv`, `references.tsv`, `target_composition.tsv` provide individual IDs and all denominators.

### AIM panel and painter

[observation] 71,552 autosomal candidates pass source divergence/call criteria;22,986 markers remain after target call>=.2 and50kb greedy thinning. Full source ALT frequencies, source/target call fractions, divergence and equal-source descriptive Fst for EVERY selected marker are saved in `aim_panel.tsv`.

| quantity | minimum | median | maximum |
| --- | --- | --- | --- |
| divergence | 0.45 | 0.639881 | 0.951509 |
| Fst_equal_source | 0.0843968 | 0.237417 | 0.746172 |
| target_call_fraction | 0.377698 | 0.76259 | 0.985612 |

[inference] Method:15 unordered ancestry-pair states;2Mb windows; diploid two-source genotype probabilities and one-allele Bernoulli pseudo-haploid probabilities; emission logs tempered by.25; row-normalized transition weights exp(-4 * minimum changed chromosome slots); uniform initial pair prior and resets between chromosomes. Forward–backward posterior chromosome fractions and Viterbi pairs are saved in `primary_paint.npz`;<5 observed SNPs produces flat emissions and excluded support. Genome-wide evaluation needs100 supported windows and20% AIM call rate. AMY support needs10 observed SNPs in supported overlapping windows. This is explicitly NOT a phased ancestry-segment assignment of H1/H2.

[inference] Known approximations: source frequencies estimated from few partially missing ancient samples; SNP independence despite LD overcounts evidence and overconfident segments;50kb thinning/.25 tempering are not LD calibration; physical-window switches are not genetic-map recombination; no genotype likelihoods, phase linkage, drift calibration or source-frequency posterior integration. A historically admixed STEPPE state and non-CHG-equivalent IRAN state bracket variation but can be mutually confounded with ancestry pairs.

### Required validation: FAILED, not waived by dominant-label accuracy

[observation] Confusion matrix: rows published source, columns inferred dominant global source; unevaluable column is excluded from accuracy. There are ZERO covered IRAN holdouts, not perfect IRAN classification.

| published_source | WHG | EHG | FARM | IRAN | STEPPE | unevaluable |
| --- | --- | --- | --- | --- | --- | --- |
| WHG | 3 | 0 | 0 | 0 | 0 | 1 |
| EHG | 0 | 4 | 0 | 0 | 0 | 0 |
| FARM | 0 | 0 | 5 | 0 | 0 | 5 |
| IRAN | 0 | 0 | 0 | 0 | 0 | 2 |
| STEPPE | 0 | 0 | 0 | 0 | 7 | 0 |

[observation] Evaluable accuracy19/19=100%, but only19/27=70.37% of held-out individuals were evaluable. Correct-source median fractions: WHG=.957974, EHG=.514091, FARM=.742684, IRAN=NA, STEPPE=.937793. Both IRAN holdouts have AIM call rates.088445 and.052945 (<.2); WHG1 and FARM5 holdouts also fail coverage. EHG median violates the prespecified.60 predominance threshold even though its largest component is EHG. The failed gates are **heldout_count, per_source_accuracy, source_dominance, stability**. Covered-subset classification cannot rescue missing-source and unstable-proportion validation.

[observation] Negative control selected WITHOUT ancestry results: highest-call-rate evaluable held-out WHG `hoe004.SG`; observed AIM call fraction=.962847, genome-wide WHG fraction=0.957974 (passes.80). PC panel contains8,937 disjoint,200kb-spaced polymorphic markers after missing/constant exclusion. R² of PCs1–3 predicted by global local-call proportions: 0.784507, 0.723711, 0.096720; mean=0.534979 (passes.50). This is descriptive in-sample association, NOT held-out ancestry truth or proof of segment accuracy. Marker coordinates/PC scores and fractions are archived.

[observation] All nine preregistered settings follow. L1 is SUM of absolute source-fraction differences, not percent agreement. Stability gate requires EVERY alternative held-out median L1<=.25 and dominant-source agreement>=.8. It fails1Mb/2 (.779726),2Mb/2 (.544155), and4Mb/2 (.250815); dominant GLOBAL labels remain unchanged but ancestry proportions do not. Target-local dominant agreement falls to43.31% for1Mb/2 despite stable global labels.

| window_bp | switch_penalty | heldout_median_L1 | heldout_dominant_agreement | target_median_L1 | target_local_median_L1 | target_local_dominant_agreement |
| --- | --- | --- | --- | --- | --- | --- |
| 2000000 | 4 | 0 | 1 | 0 | 0 | 1 |
| 1000000 | 2 | 0.779726 | 1 | 0.68585 | 1.16778 | 0.433121 |
| 1000000 | 4 | 0.121946 | 1 | 0.131028 | 0.241473 | 0.910828 |
| 1000000 | 8 | 0.162217 | 1 | 0.163523 | 0.207192 | 0.859873 |
| 2000000 | 2 | 0.544155 | 1 | 0.458775 | 0.878112 | 0.70073 |
| 2000000 | 8 | 0.233662 | 1 | 0.197767 | 0.169003 | 0.861314 |
| 4000000 | 2 | 0.250815 | 1 | 0.258772 | 0.475082 | 0.817518 |
| 4000000 | 4 | 0.0914547 | 1 | 0.0711348 | 0.0858956 | 0.912409 |
| 4000000 | 8 | 0.238158 | 1 | 0.207309 | 0.177087 | 0.839416 |

[inference] Dominant labels and PC correlations alone would misleadingly suggest success. This prototype has not demonstrated stable local ancestry proportions, much less reliable haplotype-segment backgrounds. No parameters or eligibility rules were changed to force a pass. The final-code biological rerun reproduces every finite validation quantity and sensitivity row exactly; null replaces unavailable IRAN statistics in strict JSON.

### AMY test and detection floor: deliberately NOT RUN

[observation] `amy_status.json` records STOPPED. No `amy_statistics.tsv`, local-ancestry-conditioned AMY outcome frequencies, time coefficients or P values were generated. The implemented post-gate regression/frequency/power path is unreachable in this failed biological run. **There is no valid local-conditioned detection floor to report**, because there is no validated analysis design/test.239 covered targets do not by themselves establish power after ancestry conditioning. Reporting an approximate floor based on unstable inferred fractions would violate the stop rule and invite interpretation of a broken painter.

[inference] **DOES THE AMY SIGNAL PERSIST WITHIN A FIXED LOCAL-ANCESTRY BACKGROUND? UNRESOLVED.** This local view provides no validated biological comparison and therefore DOES NOT change or corroborate the prior carriage-dominated conclusion. Power sufficiency of a valid local-background test remains unknown. Even a successful unphased dosage~local-composition regression would answer persistence conditional on estimated diploid composition, not assign either structural allele to a fixed phased chromosome ancestry. A genuine fixed-background inference requires phase/linkage validation beyond this approximation.

### Carried-forward facts independently rechecked (NOT local-painter AMY results)

[observation] Raw copy-number table + archive age metadata reproduce all supplied means on533 dated CN ancients:

| p2 | n | mean_AMY1 |
| --- | --- | --- |
| Bronze age | 46 | 6.63043 |
| CHG | 1 | 8 |
| EHG | 49 | 5.38776 |
| Early farmers | 22 | 5.45455 |
| Iron age to early modern | 333 | 7.24925 |
| Neolithic farmers | 34 | 7.14706 |
| Steppe pastoralists | 16 | 8.25 |
| WHG | 32 | 5.75 |

[observation] Independently reconstructed the prior225 matched structural diplotypes from saved cohort/dosages, refitted GLOBAL-label HC3 models, and recomputed pooled/within channel algebra. Detection floor is stated in this baseline table BEFORE interpreting its nonsignificant within effect; units are B-frequency/kyr, forward time=-BP/1000.

| cut | n | B_frequency | within_slope | HC3_SE | P | MDE80_per_kyr | pooled_slope | carriage_share |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary_le1 | 225 | 0.842222 | 0.0485754 | 0.0518258 | 0.349655 | 0.145764 | 0.0346935 | 0.946114 |
| sensitivity_le2 | 225 | 0.742222 | 0.0759689 | 0.0541781 | 0.162283 | 0.15238 | 0.0506513 | 0.942277 |

[observation] Primary previous B frequency=.842222; within=.0485754/kyr versus80%-power floor=.145764/kyr (2kyr floor=.291528), P=.349655; carriage=.946114. A<=2 within=.0759689, floor=.152380, carriage=.942277. Thus the supplied.049/.146 and~94% are primary-A<=1 quantities, not the A<=2 quantities associated with X-equivalent3.64. Saved GLMM values (artifact verification, NOT fresh GLMM fitting): A<=2 |X-equivalent|=3.637970; primary GRM retention=96.955252%. These reproduce3.64 and~97%; A<=2 clears3.16 but not5.45. Both within slopes remain below their floors; this is an underpowered baseline null, not proof of absence.

[observation] A fresh immutable released-statistics stream reproduces all canonical windows:

| locus | n | max_abs_X | above_5_45 |
| --- | --- | --- | --- |
| AMY | 2204 | 3.79496 | 0 |
| LCT | 881 | 16.5232 | 281 |
| FADS | 979 | 16.1534 | 182 |
| SLC45A2 | 1032 | 13.1496 | 221 |
| SLC24A5 | 840 | 5.99031 | 6 |

[observation] Previous matched count225/288 is verified. Raw CN+AADR mapping finds519 uniquely matched dated-coordinate ancients, versus the supplied522: **disagreement=-3**, under this explicit genetic/biological/unique-prefix mapping of533 finite-CN/age ancients. Structural mapping278 is a DIFFERENT cohort, not a replacement for the522 statement. Full mapping flags are saved. Exact canonical |X| values agree with supplied rounded values; AMY2204/3.794956889/zero above5.45 reproduces. No modern/full genome local-ancestry application was attempted.

### Reproduction and limits

[observation] Commands: `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python analysis/local_ancestry.py`; same prefix with `analysis/test_aadr_orientation.py`, `analysis/local_ancestry_verify_facts.py`, `-m unittest discover -s analysis -p test_local_ancestry.py -v`, and `analysis/local_ancestry_report.py` (APPENDS). Logs and first-run snapshots are under `data/derived/local_ancestry/`; manifests under `analysis/runs/` use `scripts/provenance.py` with all git metadata probes explicitly disabled. Source/preregistration/small-output hashes and large-input size/mtime records are provided; no large-input content hash is asserted.

[observation] Orientation reproduces low-bit-first r=-.1168, high-bit-first unflipped=-.8103, high-bit-first ALT inversion=+.8103, agreement=.7223. Six software regressions cover one-versus-two-allele weighting, pair transitions, exact-enumeration forward–backward/Viterbi, missing support/chromosome reset, conditional HC3/detection floor, prereg chronology/AIM rules/confusion counts/fail-closed outcomes. Software tests do NOT validate biological ancestry. An early coding correction ensured primary AMY coverage counts use PRIMARY rather than last sensitivity windows; final-code rerun is authoritative. An initial artifact test was explicitly skipped while biological outputs were incomplete; the completed run must pass it.

[inference] Stage1 delivers reproducible machinery and an honest failed validation, NOT a validated painter or AMY answer. A separately preregistered revision should improve reference coverage, investigate EHG/STEPPE and window/transition identifiability, and validate segment accuracy with known mosaics plus phase-linked structural calls. This run must not be tuned retroactively or scaled genome-wide.
