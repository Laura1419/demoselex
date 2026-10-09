# CN-as-SNP results

## New analysis section — incremental record, 2026-10-09

[observation] Preregistration `docs/cn-as-snp-prereg.md` was written before biological model results. The initial bounded run started at 13:34 UTC. Existing conclusions elsewhere are not rewritten.

[observation] The supervisor approved a same-cohort **ordering gate** in place of absolute reproduction; this is explicitly an amended task, not successful fulfillment of the original absolute-fidelity prerequisite. Absolute X discrepancy remains a diagnostic, and transported thresholds/extrapolated ratios cannot establish that the original scan would have called a synthetic allele.

### Increment 1 — cohort availability, before fitted results

[observation] The immutable release contains 288 ancient individuals with both categorical haplotype labels and finite dates. Exactly 225/288 (78.125%) match the released imputed-genotype bundle by unique biological-ID prefix; 63 are absent and none ambiguous. The complete mapping is `data/derived/cn_as_snp/sample_mapping.tsv`. All primary controls and, conditional on passing the validation gate, synthetic loci use the same 225 individuals. No unmatched modern reference individuals are added.

[observation] This is below the anticipated maximum n=288; no missing individuals are treated as reference homozygotes. Genome-wide GRM construction and canonical fitting are in progress. Synthetic biological trend tests have not yet been run at this increment.

## New final analysis section — 2026-10-09T13:47:56.299763+00:00

[observation] **The amended ordering gate passed.** Both canonical loci are much stronger than the released lead AMY-region real SNP in the same 225 ancients. Absolute reproduction of the paper was NOT achieved and is not asserted; the supervisor explicitly authorized the ordering gate before fitting.

| Control | Released allele (GRCh37) | Released X | Z on | X-equiv on | Z off | P ancestry on | P ancestry off |
|---|---|---|---|---|---|---|---|
| LCT | 2_136608646_G_A | 16.52320 | -5.9180 | -3.9193 | -5.9171 | 0.965 | 0.975 |
| SLC45A2 | 5_33951693_C_G | 13.14963 | -9.2970 | -6.1571 | -9.9480 | 0.305 | 0.055 |
| AMY_real_SNP | 1_104027100_G_C | 3.79496 | -0.6761 | -0.4477 | -0.6761 | 0.855 | 0.855 |

[observation] Canonical GRM-on |Z| values are 8.75x (LCT) and 13.75x (SLC45A2) the lead AMY SNP, exceeding the prespecified 1.25x separation and |Z|>=2 requirements. Both fits converged. Exact released control alleles—not regional maxima substituted at different positions—were used. ALT time-BP signs are retained: positive Z means enrichment in older samples. Signed Z differs from the released canonical X signs; the exact signed paper-slope convention is not reproduced here. ALT dosage orientation follows the supplied imputed-bundle reader, not an AADR allele substitution. Only absolute values enter the gate and extrapolation.

[observation] **Genome-wide individual-level GRM:** 1,133 polymorphic background SNPs, all 22 autosomes represented, sampled by a fixed 10,000-variant index grid then filtered for MAF>=0.05/call rate>=0.95. AMY K excludes all 75 chr1 markers, leaving 1,058 genome-wide background SNPs. LCT and SLC45A2 use the same construction rule/cohort with their own LOCO exclusion. This is not a GRM inferred from AMY and not a trajectory/group-level proxy. Marker list and actual matrices are saved. No LD pruning, paper QC recreation or paper fixed batch covariates are claimed; a sparse NoQC-bundle approximation may underrepresent structure. All control ALT dosages were integer-like (fractional counts in summary); the fitted response uses explicit ALT hardcalls, not rounded dosages.

### Synthetic loci: the substantive GRM-on versus GRM-off comparison

| Locus | Dosage n0/n1/n2 | Z on | Z off | X-equiv on | X-equiv off | P ancestry on | P ancestry off | abs(X)>3.16 on/off | abs(X)>5.45 on/off |
|---|---|---|---|---|---|---|---|---|---|
| H1^a | 164/51/10 | 4.1163 | 4.2456 | 2.7261 | 2.8117 | 0.1650 | 0.1550 | no/no | no/no |
| H5 | 145/73/7 | -2.7062 | -2.7062 | -1.7922 | -1.7922 | 0.3850 | 0.3900 | no/no | no/no |
| H7 | 186/39/0 | -1.3684 | -1.3684 | -0.9062 | -0.9062 | 0.5950 | 0.6000 | no/no | no/no |
| H9 | 223/2/0 | 1.3375 | 1.4125 † | 0.8858 | 0.9354 | 0.7771 | 0.7500 | no/no | no/no |
| H3^r | 97/85/43 | -1.4557 | -1.4606 | -0.9640 | -0.9673 | 0.0900 | 0.0800 | no/no | no/no |
| AMY1_continuous | quantitative, n=225 | -2.4804 | -2.4804 | -1.6427 | -1.6427 | 0.5700 | 0.6400 | no/no | no/no |

[observation] † H9 GRM-off **failed the primary 40-iteration fit**. The separately labelled same-model 400-iteration diagnostic converged at iteration 51 (Z=1.4125, X-equiv=0.9354, ancestry P=0.750; 199/199 diagnostic null fits succeeded). These diagnostic values do not erase the primary failure. H9 has only two alternate alleles (2/450), so numerical/inferential stability is poor. H9 GRM-on null had 165 successful and 34 failed fits out of 199 attempted; its P=0.7771 conditions on successful fits and may be biased. Every other primary converged row completed 199/199 null fits. Failed-null counts and quantiles are preserved; no failed fit is assigned Z=0.

[observation] **No structure-absorption collapse occurs in this approximation.** H1^a retains 96.96% of its GRM-off |Z| (4.1163 versus 4.2456); H5 and H7 are unchanged, H3^r retains 99.66%, and continuous AMY1 is unchanged. H9 retains approximately 94.69% relative to its longer-iteration off diagnostic. The lack of collapse with this sparse GRM is not proof that the original scan's richer GRM would behave identically. GRM-off still estimates independent overdispersion tau2; it is not a no-overdispersion binomial GLM.

### Extrapolated released-scale comparison (not measured Akbari statistics)

[observation] The prespecified number is |X_LCT,released| × |Z_h,our| / |Z_LCT,our|, with all denominators from the same 225-person cohort/model. The table uses GRM-on; GRM-off is also saved in `statistics.tsv`. The continuous row is a Gaussian-trend analogy, not a binomial allele-scan result.

| Locus | abs(Z)/LCT abs(Z) | LCT-extrapolated abs(X) | >3.16 | >5.45 | Post-hoc SLC45A2 anchor abs(X) |
|---|---|---|---|---|---|
| H1^a | 0.696 | 11.493 | yes | yes | 5.822 |
| H5 | 0.457 | 7.556 | yes | yes | 3.828 |
| H7 | 0.231 | 3.821 | yes | no | 1.935 |
| H9 | 0.226 | 3.734 | yes | no | 1.892 |
| H3^r | 0.246 | 4.064 | yes | no | 2.059 |
| AMY1_continuous | 0.419 | 6.925 | yes | yes | 3.508 |

[inference] **Uncertainty is substantial and not a numerical confidence interval.** A post-hoc alternative canonical anchor changes H1^a from 11.49 to 5.82 and H5 from 7.56 to 3.83, flipping H5's 5.45 classification. The loci do not share a validated universal multiplicative scaling factor. Small-cohort ascertainment, differing trajectories, locus-specific variances, approximation to PQL and sparse GRM all limit extrapolation. The ratio is internally comparable but not transportably calibrated; neither anchor estimates the paper's posterior pi in this cohort.

### Permutation null: ancestry-conditioned versus unrestricted

[observation] Primary shuffles permute individual dosage/phenotype values within the released **p2 ancestry component**, retaining each component's count/frequency, dates and K. Null statistics therefore need not center on zero: ancestry replacement correlated with date remains in the shuffled data. Continuous phenotype and both models are fitted anew, including variance components. An independent secondary unrestricted shuffle destroys that ancestry association. Both use seed 20261009 and 199 attempts, with plus-one two-sided empirical P. These are exchangeability nulls, not drift simulations.

| Locus/control | Ancestry P on | Ancestry P off | Unrestricted P on | Unrestricted P off | Successful unrestricted on/off |
|---|---|---|---|---|---|
| LCT | 0.9650 | 0.9750 | 0.0050 | 0.0050 | 199/199 |
| SLC45A2 | 0.3050 | 0.0550 | 0.0050 | 0.0050 | 199/199 |
| AMY_real_SNP | 0.8550 | 0.8550 | 0.4650 | 0.4650 | 199/199 |
| H1^a | 0.1650 | 0.1550 | 0.0050 | 0.0050 | 199/199 |
| H5 | 0.3850 | 0.3900 | 0.0100 | 0.0100 | 199/199 |
| H7 | 0.5950 | 0.6000 | 0.1200 | 0.1200 | 199/199 |
| H9 | 0.7771 | 0.7500 | 0.2256 | 0.1450 | 132/199 |
| H3^r | 0.0900 | 0.0800 | 0.1500 | 0.1650 | 199/199 |
| AMY1_continuous | 0.5700 | 0.6400 | 0.0250 | 0.0300 | 199/199 |
| background_chr3 | 0.8150 | 0.8200 | 0.1050 | 0.1000 | 199/199 |
| background_chr4 | 0.2300 | 0.2350 | 0.2850 | 0.2900 | 199/199 |
| background_chr6 | 0.0500 | 0.0600 | 0.0650 | 0.0650 | 199/199 |

[inference] None of the five haplotypes has ancestry-conditioned empirical P<0.05. Canonical LCT also has ancestry-conditioned P=0.965/0.975: this null asks whether the trend exceeds what survives **within these broad ancestry labels**, not whether a known selected locus has an overall temporal trend. The null must not be advertised as recovering all established selection. GRM-on/off similarity and high conditional P can coexist: a covariance random effect and conditioning permutations on categories are different operations. No strong GRM absorption here does not establish ancestry-independent selection. Individual label exchangeability is not guaranteed under kinship or correlated drift, so these empirical P values are diagnostic rather than a demonstrated kinship-aware selection calibration; no multiple-testing correction is claimed.

[observation] Full null distributions are `permutations.tsv`, `unrestricted_permutations.tsv` and the H9 diagnostic null TSV. `summary.json` and `unrestricted_null.json` give |Z| median/90th/95th/99th percentiles, completed counts and failures. The primary requested attempts all finished within the time budget; `complete_with_fit_failures` reflects numerical failures, not a timed-out null. The primary was rerun deterministically after correcting only that status label; biological statistics and nulls are independently compared to the preserved first-run snapshot.

### Optional 600-person control cohort (never used in AMY ratios)

[observation] A genotype-blind seeded sample of 600 unrelated dated ancients from Table S1 was fitted using the same marker list and LOCO construction. No synthetic AMY phenotypes are available for this cohort, so it supplies diagnostics only.

| Locus | n | Z on | Z off | X-equiv on | Released abs(X) |
|---|---|---|---|---|---|
| LCT | 600 | -10.6660 | -11.2807 | -7.0638 | 16.5232 |
| SLC45A2 | 600 | -13.3807 | -14.1052 | -8.8616 | 13.1496 |

[observation] Absolute large-|X| values are still not reproduced, and LCT versus SLC45A2 ordering differs from the full released scan. The amended gate requires separation from the real AMY SNP, not reproducing canonical ordering against each other. No full ~16k-individual fit or original scan-equivalent cohort is claimed.

### Frequencies in the primary 225 matched individuals

[observation] Ancestry frequencies (alternate copies / 2n):

| Group | n | H1^a | H5 | H7 | H9 | H3^r |
|---|---|---|---|---|---|---|
| Bronze age | 12 | 0.0417 | 0.1250 | 0.0417 | 0.0000 | 0.5000 |
| EHG | 13 | 0.2308 | 0.1538 | 0.0769 | 0.0000 | 0.3846 |
| Early farmers | 11 | 0.6818 | 0.0000 | 0.0000 | 0.0455 | 0.1818 |
| Iron age to early modern | 174 | 0.1322 | 0.2184 | 0.0948 | 0.0029 | 0.3879 |
| Neolithic farmers | 6 | 0.1667 | 0.2500 | 0.1667 | 0.0000 | 0.4167 |
| Steppe pastoralists | 3 | 0.1667 | 0.1667 | 0.1667 | 0.0000 | 0.1667 |
| WHG | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.3333 |

[observation] Date bin frequencies (alternate copies / 2n):

| Group | n | H1^a | H5 | H7 | H9 | H3^r |
|---|---|---|---|---|---|---|
| [0.0, 2000.0) | 174 | 0.1322 | 0.2184 | 0.0948 | 0.0029 | 0.3879 |
| [2000.0, 4000.0) | 10 | 0.1000 | 0.1000 | 0.0500 | 0.0000 | 0.5500 |
| [4000.0, 6000.0) | 16 | 0.1562 | 0.2188 | 0.0938 | 0.0000 | 0.3438 |
| [6000.0, 8000.0) | 14 | 0.1071 | 0.0714 | 0.0714 | 0.0000 | 0.3929 |
| [8000.0, 10000.0) | 9 | 0.6667 | 0.0000 | 0.0000 | 0.0556 | 0.1667 |
| [10000.0, inf) | 2 | 0.7500 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |

### Carried-forward facts verified from raw releases

[observation] Mean AMY1 values refer to **533 copy-number ancients**, not merely the 288 haplotype-labelled ancients or 225 matched individuals. The seven quoted values reproduce at stated precision; CHG (n=1, mean 8.00) is an additional group omitted in the carried-forward list. The exact release label is “Iron age to early modern”, not “Iron age to modern”.

| Component | Copy-number n | Mean AMY1 |
|---|---|---|
| Bronze age | 46 | 6.63043 |
| CHG | 1 | 8.00000 |
| EHG | 49 | 5.38776 |
| Early farmers | 22 | 5.45455 |
| Iron age to early modern | 333 | 7.24925 |
| Neolithic farmers | 34 | 7.14706 |
| Steppe pastoralists | 16 | 8.25000 |
| WHG | 32 | 5.75000 |

[observation] Full 288-labelled ancient haplotype component frequencies, ordered oldest to youngest by mean date:

| Component | n | Mean BP | H1^a frequency | H5 frequency |
|---|---|---|---|---|
| Early farmers | 15 | 9289.3 | 0.7667 | 0.0000 |
| WHG | 8 | 7981.4 | 0.2500 | 0.0000 |
| EHG | 14 | 6705.1 | 0.2143 | 0.1429 |
| Neolithic farmers | 6 | 5235.2 | 0.1667 | 0.2500 |
| Bronze age | 12 | 3802.7 | 0.0417 | 0.1250 |
| Steppe pastoralists | 3 | 3752.8 | 0.1667 | 0.1667 |
| Iron age to early modern | 230 | 943.7 | 0.1283 | 0.2261 |

[observation] H1^a declines in 5/6 adjacent ancestry-component steps; H5 rises from 0.0000 to 0.2261. These are component averages, not replicated within-component trajectories.

[observation] Selected modern-population checks (the broad East Asia range is evaluated on five major 1KG populations, not asserted for every tiny labelled subpopulation):

| Population | n | H5 frequency | H9 frequency |
|---|---|---|---|
| CHB | 103 | 0.3981 | 0.0194 |
| CHS | 105 | 0.4429 | 0.0286 |
| JPT | 103 | 0.5000 | 0.0097 |
| CDX | 93 | 0.4032 | 0.0215 |
| KHV | 99 | 0.4747 | 0.0303 |
| PEL | 85 | 0.4471 | 0.1176 |
| Pima | 14 | 0.7143 | 0.1071 |
| Quechua | 3 | 0.5000 | 0.1667 |

[observation] Across *all* released EA-labelled populations, including n=1–2 groups, H5 spans 0.000–0.750 and H9 spans 0.000–0.111; hence “0.40–0.50 / 0.00–0.03 in East Asia” is not literally true for every population. All per-population denominators are saved in the modern verification TSV.

[observation] Released AMY window: 2,204 variants, maximum |X|=3.794956889120, zero above 5.45, agreeing with the quoted values.

### Explicit counterfactual answer and limitations

[inference] **Under the directly fitted approximation, none of the five perfect-tag haplotype loci clears either transported cutoff** (largest |X-equivalent|=2.726 with GRM, 2.812 without). **Under the requested same-cohort LCT-anchored extrapolation, H1^a and H5 would clear 5.45** (11.493 and 7.556), while H7/H9/H3^r clear 3.16 but not 5.45. This second statement is an uncertain extrapolation, not a faithful rerun or a yes/no answer about the original Akbari scan. H1^a is 6.09x and H5 4.00x stronger in |Z| than the *released lead* AMY SNP in this cohort; the full 2,204-SNP window was not refitted to select a cohort-specific strongest tag. These contrasts are consistent with weak tagging, but are not proof that every SNP tag fails.

[inference] The original-paper question remains **unresolved**: n=225 is small and biased toward later ancients; a perfect categorical haplotype tag is an assay upper bound, not a tested variant in the real scan; ancestry-conditioned nulls do not support independent temporal change at P<0.05; sparse GRM correction does not collapse the trend but cannot settle the behavior of the paper's full relatedness correction. Increased perfect-tag signal and compatibility with ancestry-associated frequency change coexist here; weak tagging versus structure-associated change is not definitively separated for the original scan. A counterfactual threshold exceedance is not a discovery, posterior selection probability or claim that Akbari actually assayed the copy-number allele.

[observation] Reproducible entry points: `analysis/cn_as_snp.py`, `analysis/cn_as_snp_diagnostics.py`, `analysis/cn_as_snp_unrestricted.py`, and this report generator. All outputs are under allowed paths; raw releases unchanged. No git commands or staging were performed. The repository provenance harness was reused with only its git-metadata provider disabled under the explicit no-git rule; large raw genotype files have size/mtime identifiers rather than newly asserted SHA256 hashes.

## Additional validation and null-distribution summary — appended 2026-10-09

[observation] Seven new source-named tests and eight existing GLMM regressions passed. The artifact audit exactly recomputed all six saved LOCO GRMs from the same genome-wide genotype rows, verified synthetic dosages, checked X conversion and all empirical P values against their null TSVs, confirmed 4,577 primary null-attempt rows and 24 unrestricted locus/model nulls, and confirmed identical biological results across the deterministic primary rerun. No staging or git command was performed; the repository index was deliberately not inspected.

[observation] Primary ancestry-conditioned **GRM-on** absolute-Z null quantiles (GRM-off and unrestricted distributions remain fully available in the JSON/TSV artifacts):

| Locus/control | Successful / attempted | Median abs(Z) | 90th | 95th | 99th |
|---|---|---|---|---|---|
| LCT | 199/199 | 6.0619 | 6.1525 | 6.1743 | 6.1997 |
| SLC45A2 | 199/199 | 9.2233 | 9.3668 | 9.3893 | 9.4443 |
| AMY_real_SNP | 199/199 | 0.8102 | 1.0048 | 1.0570 | 1.1769 |
| H1^a | 199/199 | 3.8658 | 4.1648 | 4.2219 | 4.4172 |
| H5 | 199/199 | 2.6563 | 2.8427 | 2.8833 | 2.9317 |
| H7 | 199/199 | 1.4008 | 1.5907 | 1.6347 | 1.6958 |
| H9 | 165/199 | 1.5010 | 1.6977 | 1.7158 | 1.7319 |
| H3^r | 199/199 | 1.1849 | 1.4454 | 1.5209 | 1.7524 |
| AMY1_continuous | 199/199 | 2.5224 | 2.7774 | 2.8702 | 3.0242 |
| background_chr3 | 199/199 | 1.8817 | 2.1062 | 2.1985 | 2.3198 |
| background_chr4 | 199/199 | 0.9670 | 1.1886 | 1.2443 | 1.3314 |
| background_chr6 | 199/199 | 1.3812 | 1.5839 | 1.6232 | 1.7512 |


[observation] 174/225 (77.3%) matched individuals lie in the 0–2,000 BP bin, so the comparison is dominated by late ancient samples. The continuous phenotype is the release's per-individual integer AMY1 copy-number estimate (not an imputed mean of haplotype templates). GRM-off sets K=0; tau1 printed by the optimizer is then an inactive, nonidentified parameter, not residual relatedness.

[observation] The five major 1KG East Asian populations show H5=0.3981–0.5000 and H9=0.0097–0.0303: slight differences from the quoted 0.40 lower / 0.03 upper bounds are rounding (CHB/KHV), while the broader all-population range reported above is a genuine scope difference.
