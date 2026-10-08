# C. Akbari et al. 2026 reconstruction

> Owner task: recon-akbari. Status: IN PROGRESS.

## Paper and identifiers
- Akbari et al. 2026, Nature 654:419-428, doi:10.1038/s41586-026-10358-1
- Preprint: bioRxiv 2024.09.14.613021 / PMC11419161
- ENA project: PRJEB106907
- Code repository: https://github.com/mokar2001/PQLseqPy
- Supplementary information: 184 pages (S1-S8)

## Sample and scope

**Sample size**: 22,274 total individuals analyzed
- 6,438 modern individuals (5,935 UK Biobank SNP array, 503 1000 Genomes shotgun)
- 15,836 ancient individuals from West Eurasia
  - 5,800 previously published (3,232 capture, 2,568 shotgun)
  - 9,568 newly reported (9,487 capture, 81 shotgun)
  - 468 individuals with additional data added (181 capture, 287 shotgun)

**Geographic scope**: West Eurasia (Europe, Near East, West Asia)
- Geographical bounds: longitude between -50° and 120°, latitude north of 24°
- Defined by principal components (PC1 > 0, PC2 < -0.075, PC3 < 0.06, PC4 > -0.01)

**Temporal range**: Approximately 14,000 years BP to present
- Ancient individuals span 14,000 years
- Modern individuals provide contemporary baseline
- Uneven sampling across time with density peaks at specific periods

**Individual clustering**: Ancestry-based clustering using hierarchical clustering on top PCs
- Multiple datasets harmonized (UKBB, 1000GP, ancient individuals)
- Ancestry-matched subsets created for quality control
- Final unrelated subset: 13,936 individuals for single-variant analysis

## Pipeline

### Genotype/imputation pipeline
1. **Data sources**: Multiple ancient DNA datasets + UK Biobank + 1000 Genomes
2. **Imputation method**: IMPUTE2 with 1240k SNP panel as reference
3. **Imputation quality threshold**: IQS > 0.9 required
4. **Variant QC**: 3-step procedure with multiple filters

### Variant QC filters
**Step 1 (Provisional QC on 1240k SNPs)**
- Chi-square tests for allele frequency consistency across datasets
- INFO score > 0.6 for both shotgun and 1240k imputation
- 939,802 SNPs passed

**Step 2 (Extended QC on all imputed variants)**
- MAF > 0.002 in 1000 Genomes
- Imputation discrepancy metric: mean(|Δ|) < 0.05 across sequencing types
- Chi-square tests on ancestry-matched subsets (p > 1e-10)

**Step 3 (Batch effect correction)**
- Re-imputed UK Biobank to match ancient DNA pipeline
- GLMM with/without covariates (IQS, aDNA, SG, 1KG) for batch effect detection
- Weighted least squares outlier removal (Pearson residuals > 2)

**Final variant set**: 9,739,624 variants (18.59% of 52,382,872 imputed)
- 8,074,573 SNPs
- 1,665,051 indels

### Variants tested
All 9.7 million variants passing QC were tested for selection
- Autosomal variants only
- Includes both common and low-frequency variants (MAF > 0.002)

## The model (likelihood written out)

The GLMM uses a binomial likelihood with logit link function:

For each observation i (which may represent pooled counts from multiple individuals):

**Likelihood**:
```
Y[i,0] | p_i, n_i ~ Binomial(n_i, p_i)
```

**Linear predictor**:
```
logit(p_i) = X[i]β + u_i
```

**Random effects**:
```
u ~ N(0, V)
V = τ₁K + τ₂I
```

Where:
- Y[i,0] = count of alternative alleles
- Y[i,1] = count of reference alleles  
- n_i = Y[i,0] + Y[i,1] = total allele count
- p_i = allele frequency (latent)
- X = design matrix including time vector (calendar years BP as fixed effect)
- β = fixed effects including selection coefficient β_time
- u = random effect vector capturing population structure
- K = Genetic Relatedness Matrix (GRM)  
- I = identity matrix
- τ₁, τ₂ = variance components

**Fitting algorithm**: Penalized Quasi-Likelihood (PQL) using AI (Average Information) algorithm
- Inspired by PQLseq (Sun et al. 2019) but optimized for population genetics
- Implemented in PQLseqPy (Python)

## Selection coefficient: definition and units

The selection coefficient in this model is **β_time**, the coefficient on the time (calendar years) variable in the linear predictor:

```
β_time × time_BP
```

**Units**: Per generation change in log-odds of allele frequency per year
- Effectively represents the strength and direction of systematic allele frequency change over time
- Positive β_time indicates allele increasing in frequency over time (suggesting positive selection)
- Negative β_time indicates allele decreasing (suggesting negative selection)

**Interpretation**: This is an **inferred selection coefficient** rather than a true population genetic selection coefficient s
- Captures correlation between allele frequency and time after correcting for population structure
- Differs from classical s in discrete-generation population genetics models

## Locus-specific variance

The model includes **no locus-specific variance component** explicitly. Instead, locus-specific effects are captured through:

1. **Fixed effect β_time** - Main selection signal varying by locus
2. **Global variance components (τ₁K + τ₂I)** - Shared across all loci
3. **Residual binomial variance** - Locus-specific sampling variance

**Key point**: The model assumes that all loci share the same population structure covariance (τ₁K) but may have different time trends (β_time).

## Calibration and filters

**Calibration approach**: Uses GWAS enrichment as proxy for false discovery rate (FDR)
- Variants under directional selection should be enriched for GWAS hits
- Enrichment curve plateaus provide natural FDR control
- Specifically validated on simulations with known positives/negatives

**Significance threshold**: |Z| > 5.45 (genome-wide significance)
- Threshold determined by GWAS enrichment behavior
- More stringent than standard GWAS (|Z| > 5.45 vs ~4.9)

**Empirical validation** on simulated data shows:
- GLMM outperforms GLM with principal components across all scenarios
- Power maximal with GLMM, declines with too few or too many PCs due to stratification or collinearity

## Released data

**Public genotype/dosage data**: ENA project PRJEB106907
- Aligned sequencing data from 10,016 ancient individuals
- 367 newly reported shotgun genomes
- Capture data from 9,649 individuals (9,426 new, 223 with improved coverage)

**Summary statistics**: Available at Harvard Dataverse (doi:10.7910/DVN/7RVV9N)
- Tab-delimited compressed text file
- 9.7 million variants with GLMM statistics

**Supplementary tables**:
- Table S1: 22,274 individuals analyzed
- Table S2: 367 new shotgun genomes
- Table S3: 18,397 capture libraries
- Table S4: Selection statistics for all variants
- Table S5: 696 polygenic selection tests

**Sample dates**: Archaeological or radiocarbon dates provided as point estimates
- Broad region labels: N (North), W (West), E (East), SW (Southwest)
- Anonymized for newly reported individuals

**Individual ancestry**: Western Eurasian ancestry confirmed by:
- Principal component clustering
- Ancestry-matched quality control procedures
- Mitochondrial and nuclear contamination screening

## What is downweighted under time-GRM collinearity

When the calendar-time vector (T) is correlated with major eigenvectors of the GRM (K), the GLMM downweights the selection signal β_time to avoid confounding with population structure.

**Mathematical derivation**:

The GLMM uses the covariance structure V = τ₁K + τ₂I. When T is correlated with principal components of K:

Let K have eigen-decomposition K = UΛUᵀ, where u₁, u₂, ... are eigenvectors corresponding to population structure

When T ≈ a₁u₁ + a₂u₂ + ... (time aligned with population structure axes), the random effects component τ₁K will capture some of the temporal signal.

The estimated β_time is obtained by solving:
```
(XᵀV⁻¹X)β = XᵀV⁻¹Ỹ
```
where Ỹ is the working response.

When T is collinear with major eigenvectors of K, the effective design matrix for the fixed effects is:
```
X_eff = (I - P_K)X
```
where P_K is the projection matrix onto the column space of K.

**The collinearity effect**:
1. The time signal component that aligns with population structure gets absorbed into the random effects
2. Only the component of T orthogonal to major GRM eigenvectors contributes to β_time estimation  
3. This reduces power to detect selection when selection happened along migration/demographic axes

**Identifiability loss**: When T is fully explained by K (T ∈ col(K)), β_time becomes statistically unidentifiable - the model cannot distinguish selection from population replacement.

**Practical implication**: The GLMM is conservative - it only detects selection signals that are independent of major demographic events captured by the GRM.