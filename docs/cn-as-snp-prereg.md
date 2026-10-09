# CN-as-SNP preregistration

## New analysis section — written before results (2026-10-09)

[observation] This procedure was written before canonical-locus fitting, synthetic-locus testing, or counterfactual results. The supplied `sim/akbari_glmm.py` is explicitly an independent small-data PQL approximation, not Akbari's optimized AI-PQL implementation. No existing author's conclusions are replaced.

### Scope and approved amendment

[observation] The supervisor approved an explicit amendment before results: use the AMY-matched cohort for the primary controls and synthetic loci, with an **ordering gate** rather than requiring absolute reproduction of released X. A genotype-blind 600-person cohort is a secondary validation if time permits. Full ~16k dense fits are outside the 35-minute budget. The original absolute-reproduction requirement is still reported separately and must not be claimed satisfied by the amended gate.

[inference] Within-cohort comparisons constrain an approximation; they do not establish faithful reproduction of the paper or exchangeability between assays. Released LCT and FADS maxima are nearly tied (16.52 and 16.15), both exceeding SLC45A2 (13.15), which exceeds the AMY-window maximum (3.795). No reversed FADS/SLC45A2 ordering is assumed.

### Cohort, coding and background

[observation] Ancient H1/H2 labels and age/ancestry will be extracted from the immutable amylase v1.1 archive into `data/derived/cn_as_snp/`. Match individuals to the released imputed PLINK2 bundle by exact ID or unique biological-ID prefix before the first dot; ambiguous biological matches are excluded, and exact matches preferred. Include finite dated ancients with both haplotypes labelled. Report denominators, duplicates, matching rules, matched IDs and actual fitted n. Do not use modern samples as an unmatched baseline.

[observation] For h in H1^a, H5, H7, H9, H3^r, define y_i = 1(H1_i=h)+1(H2_i=h), n_i=2. This is a perfect biallelic tag for a categorical haplotype, not a direct biallelic representation of total copy number. Mean individual AMY1 copy number is a separate continuous quantitative phenotype, using the released per-individual estimate if available, otherwise the mean of the two labelled haplotype copy-number estimates with that approximation disclosed.

[observation] Build individual-level genome-wide GRMs from up to 10,000 regularly spaced imputed-bundle markers (systematic variant-index sampling, no phenotype or significance selection), retaining finite/polymorphic autosomal SNPs with MAF >=0.05 and call rate >=0.95. Mean-impute missing backgrounds, use `AkbariGLMM.compute_grm` VanRaden normalization, and exclude the entire tested chromosome (LOCO). For AMY, exclude chromosome 1. Never construct K from the AMY locus. Report retained marker counts and chromosomes. Dosage missing sentinels outside [0,2] become NaN. This GRM is an approximation to the paper's background/QC/relatedness construction.

### Models and statistics

[observation] Reuse `AkbariGLMM.fit`: y_i ~ Binomial(2, logit^-1(beta_0 + beta_time*time_BP_i + u_i)), u~N(0,tau1*K_LOCO + tau2*I). Fit per-locus variance components by the supplied PQL/REML routine. Genotypes must be integer hard calls for this fitter; report fractional dosages and use released hardcalls rather than silently round dosage. GRM-off uses K=0 with tau2 freely estimated (thus independent overdispersion remains); also optionally report the plain binomial GLM limit. Failed fits are failures, not significant results.

[observation] Z=beta_time/SE, X_equivalent=Z/sqrt(2.28). Positive Z means increasing toward older BP, not toward the present. Compare absolute X-equivalent to 3.16 and 5.45; these are transported reference cutoffs, not recalibrated posterior probabilities. The continuous phenotype uses Gaussian REML/GLS time slope with covariance tau1*K + sigma2*I, and K-off is ordinary independent Gaussian residuals. Its Z/X-equivalent is an analogy, not the paper's binomial scan statistic.

### Validation first and controls

[observation] Before testing any synthetic AMY phenotype, fit released SNPs at LCT/MCM6 chr2:136608646 and SLC45A2 chr5:33951693, plus the released maximum-|X| real SNP within chr1:103.5–104.5 Mb, in exactly the AMY-matched cohort, GRM on and off. Extract released X at the exact alleles; region maxima must not be mistaken for those exact SNP statistics. Additional genotype-blind background controls: first eligible SNP on each of chromosomes 3, 4 and 6 from the systematic marker grid.

[observation] Amended ordering gate passes only if both canonical SNP fits converge and both GRM-on absolute Z values exceed the AMY-region real SNP absolute Z by a factor >=1.25, with both canonical |Z| >=2. If the AMY-region real SNP is absent or fails, the gate is unevaluable and therefore fails. Do not require ordering of the two canonical SNPs against each other in this small cohort. Report exact numbers, absolute reproduction discrepancies, and the gate independently. If the gate fails or fewer than 40 AMY ancients can be matched with a usable genome-wide GRM, STOP before synthetic tests, record every synthetic statistic as not run, and state that the counterfactual is unanswered. Do not weaken the gate after seeing results.

[observation] If the gate passes, also calculate extrapolated |X| = |X_LCT,released| * |Z_AMY,our|/|Z_LCT,our| using the SAME cohort/model. This is an uncertain extrapolation, not measurement or proof that Akbari would have called the locus. Report it separately from Z/sqrt(2.28), and do not use a larger-cohort control as its denominator.

### Null and reporting

[observation] Seed=20261009. For each synthetic locus, continuous phenotype and control, shuffle whole individual phenotype values/diplotypes across individuals within ancestry component (dates and genotype K remain fixed), separately reporting unrestricted shuffling if time permits. This preserves marginal ancestry composition, dosage count distribution and date distribution; it breaks within-component phenotype–date association, not the ancestry–date relationship. Use 199 permutations per locus/model initially, incremental null rows, with empirical two-sided P=(1+# |Z_perm|>=|Z_obs|)/(1+B_success). Report requested/completed/failed permutations and null quantiles; if budget stops early report unresolved precision rather than invent P. Rare/constant phenotypes may be untestable. Conditional shuffle inference requires within-component exchangeability and does not simulate drift.

[observation] Report dosage counts and frequencies by ancestry and fixed date bins [0,2000), [2000,4000), [4000,6000), [6000,8000), [8000,10000), [10000,infinity) BP. Verify carried-forward descriptive claims only if the gate permits AMY analysis; independently verify the released AMY SNP-window count and maximum during control extraction. Results append to `docs/cn-as-snp-results.md`; preserve incremental JSON/TSV/log artifacts under permitted directories.

### Decision rule and interpretation

[inference] Strong GRM-on synthetic signal relative to canonical SNPs is consistent with weak SNP tagging masking a haplotype trend; a perfect tag is an assay upper bound conditional on this sample/model, not a guarantee of detection in the original scan. GRM-off signal falling below the reference cutoffs with GRM on is the structure-absorption case in this model; it does not identify selection versus migration. Weak signal both ways supports limited detectable change/power, not absence of selection. Failed validation means no defensible counterfactual answer. n~288 is only an anticipated maximum, to be verified.

[observation] A 35-minute session cap takes precedence over permutation completion; write partial results incrementally, never silently treat missing analyses as negative.

## Secondary diagnostic amendment — after primary fits, before diagnostic fits

[observation] Added at 13:38 UTC, after observing primary results: the primary H9 GRM-off fit failed the default 40-iteration convergence limit, with only two heterozygotes. Attempt the **same** PQL model with a 400-iteration cap as a labelled numerical diagnostic; do not erase or replace the primary failure. If it converges, evaluate 199 ancestry-stratified permutations using that same cap. This is not an alternative phenotype/model and is not a change to the passed ordering gate.

[observation] The already-approved optional secondary control cohort will be a seeded, genotype-blind random selection of 600 unrelated, dated ancient individuals from Akbari Table S1 that have exact IDs in the imputed bundle. Use the primary genome-wide marker list and the same VanRaden/LOCO rule for LCT and SLC45A2 only; never use these statistics as the denominator for the primary AMY ratio. These are diagnostic absolute-scale comparisons, not a full-cohort reproduction. Skip any unavailable secondary fit rather than change the primary cohort.

[observation] The supervisor clarified that like-for-like means the same **construction rule and cohort**, not literally identical K matrices across different chromosomes. Primary LOCO remains unchanged. No common-K sensitivity is planned within the time budget.
