
## NEW INCREMENT 1 — 2026-10-09T14:22:49.680470+00:00

[observation] Preregistration persisted before this run. Reconstructed 288 labelled ancients; exact prior cohort 225 verified. Frequency, copy-check, modern, and power/channel artifacts now saved; fitted results pending.

## NEW INCREMENT 1 — 2026-10-09T14:27:51.489938+00:00

[observation] Preregistration persisted before this run. Reconstructed 288 labelled ancients; exact prior cohort 225 verified. Frequency, copy-check, modern, and power/channel artifacts now saved; fitted results pending.

## NEW INCREMENT 2 — first completed nulls, final-code replication in progress

[observation] Detection floor before interpretation: primary 80%-power HC3 floor .145764 B-frequency/kyr (ideal independent-binomial floor .108260); observed within .048575/kyr, P=.349655, below both. Primary dosage counts n0/n1/n2=10/51/164, B=379/450=.842222. GRM-on/off Z=-4.116298/-4.245565; X-equivalent=-2.726085/-2.811694. First complete 199-draw nulls give ancestry P=.130 on/off and unrestricted P=.005 on/off. Final-code rerun reconstructs the same cohort/controls/GRMs and will preserve the full tables below.

[inference] Primary binary dosage is exactly the complement of the previous H1^a test, so this pooling cannot increase information relative to that existing binary partition. Directional pooled change is visible, but the primary within-history estimate is under the stated detection floor. GRM attenuation is small; group-centred history inference is a distinct operation.

## NEW COMPLETED ANALYSIS — 2026-10-09T14:32:26.027780+00:00

[observation] This section is appended, not a rewrite of earlier conclusions. `docs/amy-binarized-prereg.md` was written before results (local preregistration, not an external immutable registry). Primary data extraction, matching, controls, background-GRM reconstruction and both nulls completed. complete. Artifacts: `data/derived/amy_binarized/`; code: `analysis/amy_binarized.py`.

### Model choice, A/B mapping and measured-copy check

[observation] PRIMARY A = **H1^a, one haploid AMY1 copy**; B = every longer labelled haplotype. Sensitivity A<=2 = H1^a/H2A0/H2A2B2; A<=3 also adds H3^r/H3A2/H3A3B3. The initial H integer describes AMY1, not the later AMY2A/AMY2B suffixes. This cut is a **MODELING CHOICE**; the pooled B labels are not a single mutational origin. “Ancestral” is the paper framing, not a new ancestral-sequence reconstruction.

[observation] Copy checks use 4580 measured diplotypes, additive design rank 11/11, residual diploid RMSE 0.4136. Homozygote columns give measured AMY1/2 (not confounded partner means); additive coefficients estimate per-haplotype copies jointly without an intercept.

| haplotype | name_AMY1_copies | primary_le1_allele | sensitivity_le2_allele | sensitivity_le3_allele | homozygote_n | homozygote_mean_AMY1_per_hap | additive_estimated_AMY1_per_hap | diplotype_slots |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| H1^a | 1 | A | A | A | 74 | 1.15541 | 1.04314 | 949 |
| H2A0 | 2 | B | A | A | 10 | 2 | 2.02981 | 302 |
| H2A2B2 | 2 | B | A | A | 0 | nan | 2.02918 | 141 |
| H3A2 | 3 | B | B | A | 0 | nan | 3.4802 | 7 |
| H3A3B3 | 3 | B | B | A | 4 | 3.125 | 3.2526 | 81 |
| H3^r | 3 | B | B | A | 1036 | 3.00483 | 2.9875 | 4127 |
| H4A2 | 4 | B | B | B | 10 | 4.15 | 4.19832 | 228 |
| H4A2B2 | 4 | B | B | B | 7 | 4 | 4.10955 | 278 |
| H5 | 5 | B | B | B | 404 | 5.05446 | 5.05016 | 2414 |
| H7 | 7 | B | B | B | 17 | 7.26471 | 7.23346 | 529 |
| H9 | 9 | B | B | B | 1 | 9 | 9.45641 | 104 |

[inference] The integer mapping agrees approximately, NOT exactly, with measured CN. Common H3^r/H5 estimates are near 3/5; H1^a homozygote mean is 1.1554, not exactly 1; sparse H3A2 and H9 estimates are 3.4802/9.4564 versus names 3/9. H2A2B2 and H3A2 lack homozygotes. Assignment/measurement errors and linked partner composition limit the check; no labels were recoded after seeing CN.

[observation] Matched 225/288 labelled, finite-dated ancients; 63 absent, no invented baseline samples. Saved cohort order, dates and diplotypes exactly reconstructed. Genome-wide background markers 1133 across 22 autosomes, chr1-LOCO 1058. Cached chr1/2/5 matrices independently rebuilt from released background genotypes and matched to 1e-12 tolerance.

| cut | n | n0 | n1 | n2 | B_copies | alleles | B_frequency |
| --- | --- | --- | --- | --- | --- | --- | --- |
| primary_le1 | 225 | 10 | 51 | 164 | 379 | 450 | 0.842222 |
| sensitivity_le2 | 225 | 17 | 82 | 126 | 334 | 450 | 0.742222 |
| sensitivity_le3 | 225 | 96 | 96 | 33 | 162 | 450 | 0.36 |

### Detection floor — stated BEFORE interpreting trajectories

[observation] Forward time t=-BP/1000; effects are B-frequency/kyr. Approximate two-sided alpha=.05, 80% power MDE=(t_.975,217+normal_.8)*SE. HC3 charges seven ancestry baselines; independent-binomial ideal floor uses observed component B frequencies and two copies/person. All cuts have within time information 49.3137 kyr² versus pooled 1281.3338 kyr²: group-centering retains only 3.85% of time information.

| cut | B_frequency | MDE80_HC3_per_kyr | MDE80_ideal_binomial_per_kyr | MDE80_HC3_2kyr | observed_within_2kyr | date_span_kyr | MDE80_HC3_full_span |
| --- | --- | --- | --- | --- | --- | --- | --- |
| primary_le1 | 0.842222 | 0.145764 | 0.10826 | 0.291528 | 0.0971509 | 9.9565 | 1.4513 |
| sensitivity_le2 | 0.742222 | 0.15238 | 0.127613 | 0.30476 | 0.151938 | 9.9565 | 1.51717 |
| sensitivity_le3 | 0.36 | 0.11888 | 0.10774 | 0.23776 | 0.00397491 | 9.9565 | 1.18363 |

[inference] Every observed within slope is BELOW both approximate floors. Primary minimum detectable change is ~.146 B-frequency/kyr (~29.2 percentage points over 2kyr), versus observed ~.0486/kyr (~9.72 points/2kyr). The null within-history result is underpowered for the observed effect, not an exclusion of selection. Full-span MDE>1 is physically impossible under a linear frequency model and warns that a sustained full-span linear effect cannot be reliably identified with this within-component design. These Wald floors are approximate, not exact minima or kinship/drift-calibrated power; they do not say the large across-component trajectory is within-history.

### Primary trajectory: ordered ancestry components, with all counts

[observation] labelled_288 (oldest to youngest by actual mean BP):

| group | mean_BP | n | n0 | n1 | n2 | B_copies | alleles | B_frequency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Early farmers | 9289.3 | 15 | 10 | 3 | 2 | 7 | 30 | 0.233333 |
| WHG | 7981.44 | 8 | 2 | 0 | 6 | 12 | 16 | 0.75 |
| EHG | 6705.07 | 14 | 2 | 2 | 10 | 22 | 28 | 0.785714 |
| Neolithic farmers | 5235.17 | 6 | 0 | 2 | 4 | 10 | 12 | 0.833333 |
| Bronze age | 3802.67 | 12 | 0 | 1 | 11 | 23 | 24 | 0.958333 |
| Steppe pastoralists | 3752.83 | 3 | 0 | 1 | 2 | 5 | 6 | 0.833333 |
| Iron age to early modern | 943.743 | 230 | 3 | 53 | 174 | 401 | 460 | 0.871739 |

[observation] B increases in 5/6 adjacent steps; starts 0.2333 and ends 0.8717. These are donor/component means, not seven independent temporal replicates.

[observation] matched_225 (oldest to youngest by actual mean BP):

| group | mean_BP | n | n0 | n1 | n2 | B_copies | alleles | B_frequency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Early farmers | 8969.23 | 11 | 6 | 3 | 2 | 7 | 22 | 0.318182 |
| WHG | 6941.5 | 6 | 0 | 0 | 6 | 12 | 12 | 1 |
| EHG | 6383.62 | 13 | 2 | 2 | 9 | 20 | 26 | 0.769231 |
| Neolithic farmers | 5235.17 | 6 | 0 | 2 | 4 | 10 | 12 | 0.833333 |
| Bronze age | 3802.67 | 12 | 0 | 1 | 11 | 23 | 24 | 0.958333 |
| Steppe pastoralists | 3752.83 | 3 | 0 | 1 | 2 | 5 | 6 | 0.833333 |
| Iron age to early modern | 945.58 | 174 | 2 | 42 | 130 | 302 | 348 | 0.867816 |

[observation] B increases in 4/6 adjacent steps; starts 0.3182 and ends 0.8678. These are donor/component means, not seven independent temporal replicates.

[observation] sensitivity_le2 matched ancestry frequencies, ordered by mean BP; full labelled counterparts in frequencies.tsv:

| group | mean_BP | n | n0 | n1 | n2 | B_copies | alleles | B_frequency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Early farmers | 8969.23 | 11 | 8 | 1 | 2 | 5 | 22 | 0.227273 |
| WHG | 6941.5 | 6 | 1 | 4 | 1 | 6 | 12 | 0.5 |
| EHG | 6383.62 | 13 | 2 | 6 | 5 | 16 | 26 | 0.615385 |
| Neolithic farmers | 5235.17 | 6 | 0 | 2 | 4 | 10 | 12 | 0.833333 |
| Bronze age | 3802.67 | 12 | 0 | 6 | 6 | 18 | 24 | 0.75 |
| Steppe pastoralists | 3752.83 | 3 | 0 | 2 | 1 | 4 | 6 | 0.666667 |
| Iron age to early modern | 945.58 | 174 | 6 | 61 | 107 | 275 | 348 | 0.79023 |

[observation] sensitivity_le3 matched ancestry frequencies, ordered by mean BP; full labelled counterparts in frequencies.tsv:

| group | mean_BP | n | n0 | n1 | n2 | B_copies | alleles | B_frequency |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Early farmers | 8969.23 | 11 | 10 | 1 | 0 | 1 | 22 | 0.0454545 |
| WHG | 6941.5 | 6 | 5 | 1 | 0 | 1 | 12 | 0.0833333 |
| EHG | 6383.62 | 13 | 8 | 4 | 1 | 6 | 26 | 0.230769 |
| Neolithic farmers | 5235.17 | 6 | 3 | 1 | 2 | 5 | 12 | 0.416667 |
| Bronze age | 3802.67 | 12 | 6 | 6 | 0 | 6 | 24 | 0.25 |
| Steppe pastoralists | 3752.83 | 3 | 1 | 1 | 1 | 3 | 6 | 0.5 |
| Iron age to early modern | 945.58 | 174 | 63 | 82 | 29 | 140 | 348 | 0.402299 |

[observation] Matched date bins:

| group | n | n0 | n1 | n2 | B_copies | alleles | B_frequency |
| --- | --- | --- | --- | --- | --- | --- | --- |
| [10000.0, inf) | 2 | 1 | 1 | 0 | 1 | 4 | 0.25 |
| [8000.0, 10000.0) | 9 | 5 | 2 | 2 | 6 | 18 | 0.333333 |
| [6000.0, 8000.0) | 14 | 0 | 3 | 11 | 25 | 28 | 0.892857 |
| [4000.0, 6000.0) | 16 | 2 | 1 | 13 | 27 | 32 | 0.84375 |
| [2000.0, 4000.0) | 10 | 0 | 2 | 8 | 18 | 20 | 0.9 |
| [0.0, 2000.0) | 174 | 2 | 42 | 130 | 302 | 348 | 0.867816 |

### Same-cohort GLMM, reference thresholds and both anchors

[observation] Controls re-fitted from integer ALT hardcalls; prior ordering gate passed and values reproduce. PQL/REML is a small-data approximation, not the original optimized Akbari scan. GRM-off retains independent overdispersion; tau1 when K=0 is inactive. Positive Z means older BP enrichment; rising B has negative Z.

| locus | model | n | z_score | X_equivalent | released_X | converged |
| --- | --- | --- | --- | --- | --- | --- |
| LCT | GRM_on | 225 | -5.91797 | -3.91927 | 16.5232 | True |
| LCT | GRM_off | 225 | -5.91714 | -3.91872 | 16.5232 | True |
| SLC45A2 | GRM_on | 225 | -9.29696 | -6.15706 | 13.1496 | True |
| SLC45A2 | GRM_off | 225 | -9.948 | -6.58822 | 13.1496 | True |
| AMY_real_SNP | GRM_on | 225 | -0.676062 | -0.447733 | 3.79496 | True |
| AMY_real_SNP | GRM_off | 225 | -0.676062 | -0.447733 | 3.79496 | True |

| locus | model | z_score | X_equivalent | clears_3_16 | clears_5_45 | LCT_extrapolated_abs_X | SLC45A2_extrapolated_abs_X |
| --- | --- | --- | --- | --- | --- | --- | --- |
| primary_le1 | GRM_on | -4.1163 | -2.72608 | False | False | 11.4929 | 5.8221 |
| primary_le1 | GRM_off | -4.24556 | -2.81169 | False | False | 11.8554 | 5.61194 |
| sensitivity_le2 | GRM_on | -5.49321 | -3.63797 | True | False | 15.3373 | 7.76961 |
| sensitivity_le2 | GRM_off | -5.49321 | -3.63797 | True | False | 15.3394 | 7.26113 |
| sensitivity_le3 | GRM_on | -3.7141 | -2.45972 | False | False | 10.3699 | 5.25322 |
| sensitivity_le3 | GRM_off | -3.7141 | -2.45972 | False | False | 10.3714 | 4.90943 |

[observation] Primary |Z| on/off=4.1163/4.2456, retention=96.96% (3.04% attenuation), X-equivalent=2.7261/2.8117: below both 3.16 and 5.45. A<=2 gives 3.6380, above 3.16 but below 5.45; A<=3 gives 2.4597. Both alternative cuts have essentially 100% GRM retention. Neither primary nor sensitivities exhibit a structure-absorption collapse.

[inference] LCT-anchored values are 10.37–15.34 and all exceed 5.45; SLC45A2 gives 4.91–7.77, changing strict classification for A<=3. Anchor-dependent extrapolations are NOT measured released X, posterior probabilities, or evidence that Akbari would call this allele. Absolute control scale remains compressed (LCT |X-equiv|3.919 vs released16.523; SLC45A2 6.157 vs13.150).

[observation] **Primary dosage is exactly 2 minus the previous H1^a dosage**. Its primary Z magnitude and power are therefore identical to that earlier binary perfect-tag locus: pooling “all longer” has NOT created extra information relative to the H1^a-versus-rest test. It is more common than testing rare H5/H9 separately, but the earlier H1^a lane already encoded this partition. The A<=2 cut is a genuinely different binary phenotype.

### Permutation nulls (199 requested per row)

[observation] Whole dosages shuffled within p2 or unrestricted, dates/K fixed, variance components refit; plus-one two-sided P conditional on successful fits. Quantiles below are median/90th/95th/99th of absolute Z; full rows in permutations.tsv.

| null | cut | model | successful | failed | empirical_P | quantiles |
| --- | --- | --- | --- | --- | --- | --- |
| ancestry | primary_le1 | GRM_on | 199 | 0 | 0.13 | 3.8942, 4.1473, 4.2022, 4.3648 |
| ancestry | primary_le1 | GRM_off | 199 | 0 | 0.13 | 3.9848, 4.2795, 4.3749, 4.5332 |
| ancestry | sensitivity_le2 | GRM_on | 199 | 0 | 0.065 | 5.1659, 5.4390, 5.5064, 5.6133 |
| ancestry | sensitivity_le2 | GRM_off | 199 | 0 | 0.065 | 5.1804, 5.4452, 5.5064, 5.6133 |
| ancestry | sensitivity_le3 | GRM_on | 199 | 0 | 0.485 | 3.7129, 3.8952, 3.9382, 4.0548 |
| ancestry | sensitivity_le3 | GRM_off | 199 | 0 | 0.585 | 3.7412, 3.9088, 3.9714, 4.0556 |
| unrestricted | primary_le1 | GRM_on | 199 | 0 | 0.005 | 0.6678, 1.7667, 2.0235, 2.3587 |
| unrestricted | primary_le1 | GRM_off | 199 | 0 | 0.005 | 0.6745, 1.7733, 2.0208, 2.4722 |
| unrestricted | sensitivity_le2 | GRM_on | 199 | 0 | 0.005 | 0.6975, 1.7047, 1.9553, 2.5871 |
| unrestricted | sensitivity_le2 | GRM_off | 199 | 0 | 0.005 | 0.6730, 1.7317, 1.9553, 2.5871 |
| unrestricted | sensitivity_le3 | GRM_on | 199 | 0 | 0.005 | 0.5752, 1.5750, 1.8289, 2.6070 |
| unrestricted | sensitivity_le3 | GRM_off | 199 | 0 | 0.005 | 0.5779, 1.5828, 1.8795, 2.6070 |

[inference] Unrestricted null asks about total trend; the ancestry null preserves composition and asks a different question. The earlier, exactly reproduced LCT fit had ancestry P=.965 on/.975 off despite known selection: this conditioning is not a universal selection detector. Neither shuffle is a calibrated demographic/kinship-aware null. Small P unrestricted plus large P conditional is compatible with substantial carriage, not proof of selection or proof of its absence.

### Two-channel split and the effect of correction

| cut | beta_within | within_HC3_se | within_p | beta_between | between_bootstrap_se | beta_total | beta_pooled | within_pooled_contribution | carriage_absolute_share | channel_class |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| primary_le1 | 0.0485754 | 0.0518258 | 0.349655 | 0.0328241 | 0.00915167 | 0.0813995 | 0.0346935 | 0.00186948 | 0.946114 | mostly carriage/composition |
| sensitivity_le2 | 0.0759689 | 0.0541781 | 0.162283 | 0.0477276 | 0.00935233 | 0.123696 | 0.0506513 | 0.00292375 | 0.942277 | mostly carriage/composition |
| sensitivity_le3 | 0.00198745 | 0.0422673 | 0.96254 | 0.0374491 | 0.00695392 | 0.0394365 | 0.0375256 | 7.64895e-05 | 0.997962 | mostly carriage/composition |

[observation] Primary within=.048575, between=.032824, total=.081399 B-frequency/kyr. The requested total is an INDEX: actual pooled=.034694=.001869 within contribution+.032824 carriage. Carriage accounts for **94.61%** of absolute pooled contributions (94.23% A<=2, 99.80% A<=3). No cut has opposing channel signs.

[inference] **Mostly carriage/composition**, not demonstrated repeated within-history selection. A group-fixed-effect correction removes the between contribution by construction and leaves a nonsignificant within signal (primary P=.350; sensitivities .162/.963), all below detection floor. It therefore removes the statistically visible pooled trend in this operational sense, NOT 94.6% of the numerical conditional slope: the conditional coefficient increases relative to pooled while losing precision. **The genome-wide GRM random-effect correction does NOT remove the signal.** These two corrections must not be conflated. Broad period/ancestry labels are not exact histories; residual identifiability and drift calibration remain limited.

### Modern binarized replication across continents

[observation] Available populations and full cut sensitivity; n and all alternate denominators given. Every requested exact population label is present. Values for all released populations and dosage counts are in modern_frequencies.tsv.

| population | n | B_copies | alleles | B_frequency | sensitivity_le2_B | sensitivity_le3_B |
| --- | --- | --- | --- | --- | --- | --- |
| PEL | 85 | 162 | 170 | 0.952941 | 0.929412 | 0.641176 |
| Pima | 14 | 28 | 28 | 1 | 0.964286 | 0.821429 |
| Quechua | 3 | 6 | 6 | 1 | 0.833333 | 0.666667 |
| Maya | 19 | 32 | 38 | 0.842105 | 0.815789 | 0.289474 |
| Surui | 8 | 12 | 16 | 0.75 | 0.75 | 0.0625 |
| Han | 33 | 65 | 66 | 0.984848 | 0.969697 | 0.530303 |
| CHB | 103 | 196 | 206 | 0.951456 | 0.936893 | 0.504854 |
| CHS | 105 | 194 | 210 | 0.92381 | 0.914286 | 0.533333 |
| JPT | 103 | 199 | 206 | 0.966019 | 0.961165 | 0.572816 |
| CEU | 99 | 168 | 198 | 0.848485 | 0.727273 | 0.333333 |
| TSI | 107 | 193 | 214 | 0.901869 | 0.831776 | 0.364486 |
| IBS | 107 | 192 | 214 | 0.897196 | 0.841121 | 0.341121 |
| GBR | 91 | 155 | 182 | 0.851648 | 0.747253 | 0.296703 |
| YRI | 108 | 197 | 216 | 0.912037 | 0.842593 | 0.305556 |
| GWD | 113 | 196 | 226 | 0.867257 | 0.80531 | 0.256637 |
| ESN | 99 | 180 | 198 | 0.909091 | 0.893939 | 0.30303 |
| Mbuti | 14 | 28 | 28 | 1 | 0.857143 | 0.285714 |

[observation] Primary PEL .9529, Pima/Quechua 1.0, Han .9848 and CHB/CHS/JPT .9238–.9660 exceed Surui .7500 and most Europe .8485–.9019. However Maya .8421 is lower than several European/African samples, YRI .9120/ESN .9091/GWD .8673 are already high, and **Mbuti is 1.0** (28/28 alleles).

[inference] Descriptive high-B in several agricultural populations on separate continents is present, but **not a specific or universal agricultural replication under the primary cut**: the hunter-gatherer comparator Mbuti saturates too. The higher-copy cuts sharpen some contrasts (see table) without establishing independent selection origins. Shared labels, ancestry, small Pima/Quechua denominators and postulated starch categories preclude a causal diet claim; B-commonness may erase high-copy functional distinctions.

### Cut sensitivity and explicit first-question verdict

[inference] **YES, a pooled directional binary AMY trend is visible** under the prespecified unrestricted-P rule, consistently across cuts. **NO, the sparse LOCO GRM does not remove it.** **YES, ancestry fixed-effect/within-history conditioning removes statistical detectability, leaving mostly carriage and an underpowered within estimate.** Primary directly fitted X-equivalent does not clear either transported threshold; the <=2 sensitivity clears 3.16 only. Thus direction and correction conclusions are robust, but cutoff classification and modern dietary specificity depend on the MODELING CHOICE. This is not a faithful original-paper counterfactual or new selection discovery.

## SECOND QUESTION — separate descriptive released-|X| distribution

[observation] Independently streamed 9,739,624 finite released variants, 0 missing X; no centering/FILTER restriction. |X| quantiles:

| quantile | abs_X |
| --- | --- |
| 0 | 1.80348e-08 |
| 0.25 | 0.461365 |
| 0.5 | 0.981059 |
| 0.75 | 1.696 |
| 0.9 | 2.4705 |
| 0.95 | 2.99238 |
| 0.99 | 4.26717 |
| 0.999 | 7.24272 |
| 1 | 16.6433 |

[observation] Fraction strictly exceeding each canonical exact window maximum (ties excluded), alongside operator-rounded thresholds:

| locus | abs_X | n_exceed | fraction_exceed |
| --- | --- | --- | --- |
| LCT | 16.5232 | 2 | 2.05347e-07 |
| FADS | 16.1534 | 4 | 4.10693e-07 |
| SLC45A2 | 13.1496 | 583 | 5.98586e-05 |
| SLC24A5 | 5.99031 | 22110 | 0.00227011 |
| AMY | 3.79496 | 167199 | 0.0171669 |
| rounded_LCT | 16.52 | 3 | 3.0802e-07 |
| rounded_FADS | 16.15 | 5 | 5.13367e-07 |
| rounded_SLC45A2 | 13.15 | 583 | 5.98586e-05 |
| rounded_SLC24A5 | 5.99 | 22115 | 0.00227062 |
| rounded_AMY | 3.795 | 167193 | 0.0171663 |

[observation] LCT16.5232 versus AMY3.79496 is a **12.7282-X gap, 4.3540-fold magnitude**; only 2 variants exceed canonical LCT versus 167,199 (1.71669%) exceeding AMY. SLC24A5 is already near the strict threshold at5.9903. Genome-wide maximum is 16.643323946 at 5_131664516_TTATATA_T, so the operator’s claim that LCT is the strongest signal is not literally true at the released single-variant level: this chr5 variant is larger than canonical LCT16.52. Top20 IDs are in distribution.json.

| category | variants | variant_fraction | any_blocks | max_class_blocks |
| --- | --- | --- | --- | --- |
| 3.16<=|X|<=5.45 | 356335 | 0.0365861 | 10007 | 9285 |
| |X|>5.45 | 31939 | 0.00327928 | 723 | 723 |

[observation] 250kb blocks use (chr,floor(POS/250000)), 10672 represented. Any-band and any-strong sets overlap in 722 blocks; mutually exclusive block-MAX categories are 9285 band-only versus 723 strong. These are not independent LD units.

[inference] The strict threshold keeps **31,939 variants in 723 blocks**, including all listed canonical windows except AMY; it discards **356,335 band variants**, including AMY’s best SNP and variants in **9,285 additional band-only blocks**. The band contains ~11.16 times as many variants as the retained strong set. There is a broad prioritization domain below the strict gate, not just the handful of spectacular LCT/FADS signals. **False-positive caveat: the discarded band is NOT all real selection**; released FDR/inflation, demography, LD and correlated tests make these follow-up candidates, not discoveries. Magnitude proximity is descriptive, not shared mechanism or causal evidence.

### Carried-forward verification and disagreements

[observation] Measured CN means/slopes rechecked on the 533 dated CN ancients, NOT the 225 matched or 288 labelled cohort:

| component | n | mean_AMY1 | beta | p |
| --- | --- | --- | --- | --- |
| Bronze age | 46 | 6.63043 | 0.433137 | 0.602579 |
| CHG | 1 | 8 | nan | nan |
| EHG | 49 | 5.38776 | 0.148729 | 0.434315 |
| Early farmers | 22 | 5.45455 | 0.407479 | 0.192887 |
| Iron age to early modern | 333 | 7.24925 | -0.0572451 | 0.944731 |
| Neolithic farmers | 34 | 7.14706 | -0.338942 | 0.665146 |
| Steppe pastoralists | 16 | 8.25 | 0.0613734 | 0.955355 |
| WHG | 32 | 5.75 | 0.0427397 | 0.855244 |

[observation] All identified CN within-component slopes are nonsignificant, not proved zero; CHG n=1 is unidentifiable. H1^a full labelled component frequency .7667→.1283, declining 5/6 steps, and H5 .0000→.2261 reproduce; verification table follows.

| component | n | mean_BP | H1a_frequency | H5_frequency |
| --- | --- | --- | --- | --- |
| Early farmers | 15 | 9289.3 | 0.766667 | 0 |
| WHG | 8 | 7981.44 | 0.25 | 0 |
| EHG | 14 | 6705.07 | 0.214286 | 0.142857 |
| Neolithic farmers | 6 | 5235.17 | 0.166667 | 0.25 |
| Bronze age | 12 | 3802.67 | 0.0416667 | 0.125 |
| Steppe pastoralists | 3 | 3752.83 | 0.166667 | 0.166667 |
| Iron age to early modern | 230 | 943.743 | 0.128261 | 0.226087 |

| population | n | H5_frequency | H9_frequency |
| --- | --- | --- | --- |
| PEL | 85 | 0.447059 | 0.117647 |
| Pima | 14 | 0.714286 | 0.107143 |
| Quechua | 3 | 0.5 | 0.166667 |
| Maya | 19 | 0.184211 | 0.0526316 |
| Surui | 8 | 0.0625 | 0 |
| Han | 33 | 0.454545 | 0 |
| CHB | 103 | 0.398058 | 0.0194175 |
| CHS | 105 | 0.442857 | 0.0285714 |
| JPT | 103 | 0.5 | 0.00970874 |
| CEU | 99 | 0.222222 | 0 |
| TSI | 107 | 0.163551 | 0 |
| IBS | 107 | 0.186916 | 0.0046729 |
| GBR | 91 | 0.153846 | 0.00549451 |
| YRI | 108 | 0.217593 | 0 |
| GWD | 113 | 0.181416 | 0 |
| ESN | 99 | 0.227273 | 0 |
| Mbuti | 14 | 0.25 | 0 |

[observation] Disagreements with supplied modern approximations: PEL H5=.44706, NOT .46 (rounds .45); Pima=.71429, Quechua=.50. CHB H5=.39806 and KHV H9=.03030 marginally differ from .40/.03 at precision; small EA groups are not covered by those broad ranges. PEL/Pima/Quechua H9=.11765/.10714/.16667 agree with .11–.17 after rounding. AMY released window remains 2,204 SNPs, max3.794956889, zero >5.45. Full released band356,335 is six more than prior two-channel scored moderate356,329: different sampled panel, not an inconsistency in the released counts. Full represented blocks10,672 versus prior scored panel10,634 likewise differ in ascertainment. Canonical rounded values agree at stated two decimals; exact and rounded exceedance counts differ, as explicitly tabulated.

### Validation, provenance and limits

[observation] Source-named tests audit dosage/complement identity, all count denominators, three cut definitions, within/pooled channel algebra, control reproduction, preregistration-before-output timestamp, null rows/P values, and genome-wide block totals. Commands and validation log are archived. Provenance manifests under analysis/runs/ use scripts/provenance.py with ONLY git metadata probes disabled. Raw inputs untouched; no git commands or staging performed. Large genotype/stat files recorded by path/size/mtime, not falsely asserted content hashes.

[inference] Sparse GRM, n=225 late-biased sampling (174/225 below2kyr), label coarseness, copy/genotype uncertainty, approximate Wald power, uncalibrated shuffle/drift nulls and anchor-dependent scale limit interpretation. Directional pooled frequency change is measured; independent within-history selection and what the full original Akbari GRM would do remain unresolved.

## NEW VALIDATION ADDENDUM — 2026-10-09T14:33Z

[observation] Seven new `analysis/test_amy_binarized.py` tests and eight existing GLMM pytest regressions pass. The final-code biological run reproduces the first completed run exactly for controls, all fitted allele statistics, channels, counts, and all **2,388 permutation rows** (12 rows × 199 successful draws, zero failures). Genome-wide streaming independently rechecks every canonical window maximum, including AMY 2,204/3.794956889/zero strict hits. All three final-code cached LOCO matrices were independently rebuilt from raw genotype markers. An interim audit run during ongoing permutations correctly failed its completeness assertion; that was partial-run state, not a final regression failure. An initial extraction/reporting run encountered a duplicate `n` dictionary keyword; corrected before fitted results and both completed deterministic runs.

[observation] Reproduction: `PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 .venv/bin/python analysis/amy_binarized.py --mode biology --seconds 850`; same entry point `--mode distribution`; `.venv/bin/python analysis/amy_binarized_report.py` (APPENDS a new section); `.venv/bin/python -m unittest discover -s analysis -p test_amy_binarized.py -v`; `.venv/bin/python -m pytest sim/test_akbari_glmm.py -q -p no:cacheprovider`. Preserve `PYTHONDONTWRITEBYTECODE=1` for imports outside allowed write roots. Main summaries/TSVs, logs, first-run snapshots and provenance hashes are saved. No raw modifications, git commands, staging or commits were performed; pre-existing index contents were not inspected under the no-git directive.
