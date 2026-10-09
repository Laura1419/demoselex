# Two-channel results — NEW APPENDED ANALYSIS, 2026-10-09

[observation] `docs/two-channel-prereg.md` was written BEFORE results and before running `analysis/two_channel.py`. This document is a new analysis, not a replacement of previous conclusions. Seed 20261009; correlated-null Monte Carlo 4096 draws/scenario; AMY member bootstrap seed 20261010, 1000 replicates. All artifacts live under `data/derived/two_channel/`; complete SNP table is `data/derived/two_channel.tsv`.

## AMY motivating positive example — first incremental result

[observation] Independently extracted the Bolognini/Sudmant/Garrison analysis archive to `analysis/mixture_null_scratch/` and joined ancient sample CN to numeric archive age. N=533, eight measured labels including CHG n=1 (CN=8). The seven supplied means reproduce: WHG 5.7500 (n32), Early farmers 5.4545 (n22), EHG 5.3878 (n49), Neolithic farmers 7.1471 (n34), Steppe pastoralists 8.2500 (n16), Bronze age 6.6304 (n46), Iron age to early modern 7.2492 (n333). Every identified component slope is nonsignificant (HC3 P=.193–.955); this is NOT proof they are zero.

[observation] Both channels, copies/kyr: **beta_within=.151399**, HC3 SE=.119441, P=.205519; **beta_between=.207252**, member-stratified bootstrap SE=.038117; requested additive **beta_total=.358650**. The numerical statement “beta_within~0” means indistinguishable from zero, NOT an exactly zero estimate; our estimate is .151399. Both remain in the ranking (inclusive AMY score 4.37791); carriage is not subtracted or demoted.

[inference] The additive total is an index, not the observed pooled slope: pooled observed=.218218, within pooled contribution=.010966, and .010966+.207252=.218218. Group-centering leaves little pooled time information; conflating the .151399 conditional slope with its .010966 pooled contribution would exaggerate the observed total. Neither constituent is discarded. Between-channel change is an important, measured part of this motivating example.

[observation] Ordered by actual MEAN date, the seven >=8-person donors are Early farmers → WHG → EHG → Neolithic farmers → Steppe → Bronze age → Iron/early modern. This differs from simply arranging the labels WHG first: Early farmers have the earlier mean date (8414 versus 7694 BP). Four of six successive steps increase. Spearman rho=.714286; relatedness-discounted signed concordance C=.582537. N_eff=2.00642, not seven independent donors. Genetic-covariance drift P=.107152; half shared-history P=.195753; fully cumulative shared-history drift P=.353429; prespecified conservative ordering P=.353429. These are actual simulated correlated drift nulls, not an independence permutation or an Akbari-derived P.

[inference] The upward tendency is real descriptively, but **does not reject this preregistered correlated/nested drift null**. This is not a dismissal of between-channel selection: significant robust ordered carriage WOULD count as candidate selection evidence. AMY is retained as a motivating phenotype/CN positive example with substantial carriage, shared haplotypes, and uncertain source mechanism, not claimed as a newly significant ordering detection. An insufficient ordering result cannot exclude selection in forming/unsampled source populations or CN-specific selection. West Eurasian donors are partly nested: steppe inherits CHG/EHG, and period labels are not independent lineages. The measured trajectory-proxy covers 50% of steppe members and 77–98% of other nonsingleton groups; it is not a fitted pedigree. Adding cumulative shared drift expressly confronts neutral inheritance, rather than counting each incoming population as a fresh independent replicate.

[observation] Global replication context: PEL mean AMY1=9.51765 (n85), Pima=10.28571 (n14), Surui=5.25000 (n8), agreeing with quoted 9.52/10.29/5.25 at rounding. H5 slot frequency East Asia=.43602 (1438 slots; quoted .44), WEA=.18642 (1738), five-population Europe=.17097 (1006), PEL=.44706 (170), Pima=.71429 (28). Explicit discrepancies: **PEL .44706, not .46**; **Europe .17097, not .18** for CEU/FIN/GBR/IBS/TSI. WEA includes Western Asia and is not identical to Europe. Pima rounds .71 and East Asia .44. Shared H5 means carriage is a live alternative; a shared structural label alone does not establish identity of mutational origin or causal ancestry. Modern high CN is replication context, not an independent ancient selection test.

## Genome-wide / Akbari audit

[observation] In progress; appended once the registered full run and independent released-statistics stream finish. No genome-wide biological inference is made from partial output.

## Carried-forward limits

[observation] Agreement with Akbari is NOT independent validation: same imputed genotypes. Earlier top-N bootstrap retention was 43.7%; no new rank-stability estimate is substituted here. A calibrated LD-honest block drift-null was NOT achievable with the existing simulation and remains unavailable.

[inference] The correlated Gaussian null is explicit and conditional, but neither a calibrated full-demography null nor a CN evolutionary-process model. Controls are sampled loci, not established neutral loci. False positives are acceptable for phenotype prioritisation; actual control fractions, not significance vocabulary, determine the practical follow-up burden. 250kb blocks and matching control-only block fractions will be reported, not variant counts alone.

## Completed genome-wide run — appended 2026-10-09T02:56Z

[observation] All eight requested merged panels were processed, yielding **538,262 unique scored variants**: 31,939 Akbari-strong, 356,329 moderate, and **149,994 controls**, not exactly 150,000. The quoted 150k is six more than the scored control set; no additional variants were invented. Exact source-locus coverage is checked against all eight `priority_loci.tsv` files. Control coverage pass fraction is 149,994/149,994=100%. The uncentered released X is used only for panel labels, not ordering or slope fitting. |X|>5.45 is the supplied pi>99% criterion; 3.16 is the approximate 50%-FDR boundary, not proof of neutrality.

[inference] SNP units are ALT-frequency/kyr. beta_within reuses the measured, within-deme-centred estimates; beta_between reuses the measured-member ancestry-turnover slope and is independently regression-checked on cached SNP genotypes. Donor means for ordering are freshly recomputed. beta_total is their requested sum, not the pooled regression. No term is zeroed, and opposite signs do not cancel ranking evidence. Between-channel background scale is control RMS **.00623638/kyr**; “large between” is |beta_between|>=.01247276/kyr. “Large within” is |beta_within|/existing corrected_SE>=2. This comparison includes a TRAIN-fitted shared within temporal-variance law carried from the earlier scan; the between RMS is estimated on all sampled controls and is a descriptive magnitude yardstick, not held-out validation.

[observation] The genome-background relatedness proxy mapped 2,210/2,265 panel individuals (97.57%) by biological ID. It averages old-trajectory covariance into rebuilt demes; it is neither a new independent pedigree nor fractional donor ancestry. Null draws use measured donor counts/dates and the same discounted pair statistic at each locus. Each locus reports `beta_within`, `beta_between`, `beta_total`, signed `ordering_statistic`, Spearman `ordering_rho`, successive-step fraction, three scenario P values, conservative `ordering_p`, donor count/N_eff, both magnitude scores, inclusive score and global rank. The mapping and covariance matrices are archived alongside the table. Whole-scan channel/control results follow; the matching control fraction is reported for EACH criterion, not a single control number reused as an omnibus calibration claim.

| Criterion | All-variant count / fraction | Matching control count / pass fraction (n=149,994) |
|---|---:|---:|
| Carriage class: large between, small within | 75,876 / 14.0965% | 9,303 / **6.2022%** |
| Within-deme class: large within, small between | 60,186 / 11.1815% | 5,182 / **3.4548%** |
| Both large | 29,436 / 5.4687% | 469 / **0.3127%** |
| Neither large | 372,764 / 69.2533% | 135,040 / **90.0303%** |
| Opposing slope signs (not necessarily large) | 133,532 / 24.8080% | 67,658 / **45.1071%** |
| Robust ordered series, P<=.05 | 93 / .01728% | 0 / **0%** |
| Robust ordering AND large between (including both-large) | 93 / .01728% | 0 / **0%** |
| Robust ordering AND carriage-class (small within) | 0 / 0% | 0 / **0%** |
| Inclusive top decile | 53,827 / 10.0001% | 1,242 / **.8280%** |
| Top decile AND ordered carriage-class | 0 / 0% | 0 / **0%** |

[inference] The 75,876 carriage-class loci are a substantial retained prioritisation channel (matching control pass 6.2022%), not excluded nuisance. They are NOT 75,876 discoveries or a percentage of selection explained. The 93 robust ordered large-between candidates all also have large within change and Akbari-strong labels (matching ordering-control fraction 0/149,994=0%); **large between with significant correlated-null ordering is candidate selection evidence here**, even when a within signal also exists. This scan found no robust ordered large-between/small-within class (matching control fraction 0%). Thus it does NOT establish an extra low-Akbari selection class from ordered carriage under this null; it also does not make carriage vanish. Strong-only ordering could reflect assay/LD/model properties and score-label reuse, not independent validation.

[inference] Inclusive ranking is `log(1+z_within²)+log(1+z_between²)+max(0,-log10(P_order)-log10(20))`. The two channel magnitudes are additive and equally capable of raising rank; ordering adds a bonus. Akbari agreement is descriptive only. Large carriage without robust ordering, within change or separate functional evidence remains a follow-up hypothesis, not selection evidence by itself. This is the explicit falsifier for claiming this measured ordering establishes more than the registered drift null; it does not falsify selection in an unobserved ancestor.

### LD-aware reporting: 250kb block representatives

[observation] There are **10,634 represented 250kb blocks**, ranked by maximum inclusive-scoring representative, with **626 control-only blocks** (all sampled variants in that block have |X|<3.16). Control-block coverage pass fraction 626/626=100%; represented control loci have coverage 100%. These are block representatives, not an LD-honest block P-value correction. Mixed blocks are not designated neutral merely because they contain a control SNP. `data/derived/two_channel/blocks.tsv` reports block ranks and both channels/ordering for each representative.

| Representative criterion | Block count / fraction | Control-only block count / fraction (n=626) |
|---|---:|---:|
| Carriage class | 3,520 / 33.1014% | 146 / **23.3227%** |
| Within class | 2,033 / 19.1179% | 81 / **12.9393%** |
| Both large | 992 / 9.3286% | 11 / **1.7572%** |
| Neither large | 4,089 / 38.4521% | 388 / **61.9808%** |
| Opposing signs | 3,149 / 29.6126% | 274 / **43.7700%** |
| Robust ordering | 4 / .03762% | 0 / **0%** |
| Robust ordering AND large between | 4 / .03762% | 0 / **0%** |
| Ordered carriage-class | 0 / 0% | 0 / **0%** |
| Inclusive block top decile | 1,064 / 10.0056% | 11 / **1.7572%** |
| Top decile AND ordered carriage-class | 0 / 0% | 0 / **0%** |

[inference] The variant-level inclusive top-decile follow-up proxy is 1,242 sampled control hits (control pass .8280%), **2.307% of the 53,827 follow-ups**. At block level, 11 control-only hits (control pass 1.7572%) are **1.034% of 1,064 follow-ups**. These are accepted practical false-positive-cost proxies, NOT estimated FDR: controls may contain genuine selection, mixed candidate/control blocks are excluded from the control-only denominator, panel ascertainment differs, and maximum selection depends on variants/block. The much higher carriage-class control-block fraction (23.3227%, vs variant 6.2022%) illustrates the cost of maximum-over-many-variants screens. Zero observed ordering-control hits does not prove zero false positives or calibrated discovery error; locus LD invalidates an independent-binomial interpretation.

### Independent released Akbari statistics audit

[observation] Streamed the immutable full statistics release with `analysis/load_akbari.py`, rather than the prioritisation subset. chr1:103.5–104.5Mb contains **2,204 SNPs**, maximum **|X|=3.794956889**, locus `1_104027100_G_C`, and **zero** |X|>5.45. This reproduces the operator's 2,204/3.795 result at rounding: released SNP statistics do NOT detect AMY at that threshold. For this audit the corresponding genome-wide ordering-control pass remains 0/149,994=0%; the local X audit itself is not a new genome-wide discovery test.

[observation] Canonical-window maxima (exact output and all bounds in `akbari_audit.json`): LCT 16.523199594 (881 SNPs), FADS 16.153410398 (979), SLC45A2 13.149626230 (1,032), SLC24A5 5.990305093 (840). Disagreements with the approximate quoted numbers are explicit: **FADS 16.15341 versus 16.1** (rounds 16.2 at one decimal), **SLC45A2 13.14963 versus 13.2** (rounds 13.1 at one decimal). LCT agrees with 16.5 and SLC24A5 with 6.0 at one decimal. X was never centred. This contrasts an undetected AMY SNP window with CN evidence; it does not establish that Akbari's SNP statistic has directly tested AMY1 copy number. Matching ordering-control pass fraction is 0% as above; canonical X maxima are released-statistic observations, not independent validation of this ranking.

### What this does / does not establish

[observation] The requested two channels and ordered-series test exist for every scored locus, and both enter the ranking. AMY has large measured between-channel CN change despite no released SNP threshold detection. Its mean-dated donor progression is upward but does not survive the explicit correlated/nested drift test (P=.35343).

[inference] This establishes an inclusive **candidate prioritiser**, not a causal partition or a new validated selection scan. The 93 ordered-large-between loci (control pass 0%) are candidates conditional on this Gaussian covariance model; no calibrated LD-honest block drift-null exists here. None is claimed as independent Akbari validation. No SNP-to-CN variance conversion, selection coefficient estimate, independent mutation count, or universal rejection of AMY selection is made. Carried-forward 43.7% top-N bootstrap retention warns that precise ranks remain unstable. Independent ancestry/CN haplotypes, phenotype follow-up, and LD/demography calibration remain necessary.

## Final validation and descriptive block cross-check — appended 2026-10-09T02:59Z

[observation] `analysis/validate_two_channel.py` passes: exact 538,262-locus source-set equality, unique IDs, finite per-locus channels/ordering, beta_total identity, inclusive score identity, descending rank, unique 10,634 block representatives, AMY pooled scaling, and 2,204/3.794956889 audit. Five new software tests pass, including baseline invariance, actual cached-versus-recomputed SNP slopes/composition, correlated-null determinism, orientation and nesting discounts. Seven existing prioritisation tests also pass. Neither regression suite fits a scientific null or validates LD calibration.

[observation] To distinguish maximum-score block representatives from an “any ordered locus” block screen, the 93 ordered variants span **eight** 250kb blocks (control-only pass 0/626=0%; corresponding variant-control pass 0/149,994=0%). Only four have an ordered maximum-score representative, as reported above. `ordered_blocks.tsv` lists the best-ranked ordered locus per such block, retaining both channels. These are candidate representatives, not eight independent discoveries:

| Locus | beta_within | beta_between | beta_total | C | correlated-null P |
|---|---:|---:|---:|---:|---:|
| 6_33161660_T_TG | -.054534 | -.047254 | -.101788 | -.731552 | .046131 |
| 6_31945646_C_T | -.050932 | -.054849 | -.105781 | -.739559 | .039053 |
| 5_131436896_C_T | .049763 | .047376 | .097139 | .728796 | .048328 |
| 11_88549245_A_G | -.035761 | -.050172 | -.085934 | -.744482 | .034904 |
| 11_88485142_G_A | .033588 | .049888 | .083475 | .744224 | .035148 |
| 4_38832914_C_CAA | -.040375 | -.026932 | -.067308 | -.730593 | .047352 |
| 6_28092993_T_TTTTTG | .030838 | .019134 | .049972 | .733007 | .044911 |
| 6_28825286_C_CT | .024084 | .021712 | .045796 | .761459 | .023920 |

[inference] The merged scored panel contains indels as well as SNPs (visible in allele IDs above); these are variant-frequency channels, not copy-number estimates. Adjacent 250kb bins may remain LD-correlated (notably chr6 and paired chr11 bins), so eight blocks are NOT eight independent donor/selection replications. Control-only block pass fraction is 0%, not a calibrated block FDR.

[observation] Inclusive top-decile ranking retains **16,822 carriage-class** variants (matching control intersection 655/149,994=.4367%), 7,959 within-class (136/149,994=.09067%), and 29,046 both-large (451/149,994=.30068%). This directly demonstrates that the BETWEEN channel contributes to ranking rather than being used to cancel the WITHIN channel. Source-defined control class alone has 1,242 top-decile hits (.8280% control pass). Exact counts/intersections are persisted in `validation.json`.

[inference] `validation.json` additionally gives nominal independent-binomial intervals for zero ordering-control counts: variant [0,.00002459] and block [0,.00587545]. They are arithmetic sampling summaries ONLY: LD and ascertainment invalidate treating them as discovery-error guarantees. The calibrated LD-honest block null remains unavailable.

[observation] Reproduction commands (run with `OPENBLAS_NUM_THREADS=1` to bound resources): `.venv/bin/python analysis/two_channel.py --mode amy`; `--mode audit`; `--mode genomewide`; `.venv/bin/python analysis/validate_two_channel.py`; `.venv/bin/python -m unittest discover -s analysis -p test_two_channel.py -v`; the existing `test_prioritization_score.py` suite. The full SNP and AMY modes were run twice with the recorded seeds; scientific estimates were unchanged. Provenance uses `scripts/provenance.py` with all git probes disabled, hashes source/small input/output files and records large immutable inputs by path/size. No raw file was edited and no git command, staging or commit was performed. Preregistration persisted 02:49:12Z, before first two-channel AMY output at approximately 02:51Z. Results are a local pre-registration, not an external immutable registry.

[observation] Final reproducibility cross-check: output SHA256 hashes match across the two AMY runs for amy_summary.json and across both full SNP runs for all eight per-panel channel tables. No statistical result changed on repetition. Final source/docs syntax, whitespace and conflict-marker checks pass without git. Tests were re-run once after the first validation manifest; the final validation manifest captures the current logs rather than rewriting the earlier record. Complete implementation finished well within the approximately35-minute bound.
