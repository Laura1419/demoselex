# B. Literature map

> Owner task: recon-lit. Status: COMPLETED, upgraded by citation-verification pass.
> Survey date: 2026-10-08. Citations were verified against primary sources via web
> search during the 2026-10-08 passes; any citation we could not fully verify is
> marked "(unverified)".

This document surveys the closest existing method families to the proposed approach —
inferring temporal population trajectories from ancient DNA and testing *replicated
directional allele-frequency change* across them. Each family is scored against the
seven dimensions required by the task:

1. **Unit of analysis** — the fundamental element being tested
2. **Structure representation** — how population structure/demography enters the model
3. **Nuisance** — which aspects of structure are treated as confounding noise
4. **Selection evidence** — the statistical criterion for calling selection
5. **Per-population selection?** — can selection differ among populations
6. **Repeated directionality?** — can the method detect the same direction of change
   across multiple distinct histories as evidence
7. **Adaptable here?** — how (if at all) the method transfers to our problem

A note on naming: the task list mentions "CLES/CLUES". We could not find any distinct
published population-genetics method named **CLES** (searched 2026-10-08); the only
methods matching are **CLUES** (Stern et al. 2019) and its extension **CLUES2**
(Vaughn & Nielsen 2024). We therefore treat "CLES/CLUES" as one family and split out
"ancestry-stratified CLUES" as an emerging practice (row 7). If "CLES" refers to a
tool used inside the Akbari et al. (2026) supplement, that should be confirmed
against the supplement — recorded as an open verification item.

## Methods summary table

| Method family | Unit of analysis | Structure representation | Nuisance | Selection evidence | Per-population selection? | Repeated directionality? | Adaptable here? |
|---------------|------------------|--------------------------|----------|-------------------|---------------------------|--------------------------|-----------------|
| **1. Ancient-DNA allele-frequency time-series selection** (Mathieson et al. 2015; Foll et al. 2015 WFABC; Akbari et al. 2026) | Locus in a dated sample series, pooled or grouped by label | Fixed group labels (period/culture/region) or genetic PCs/GRM as covariates in a regression/GLMM | Population structure, migration, turnover (removed via covariates); ascertainment, damage, coverage | Consistent trend of allele frequency (or dosage) on time; ABC posterior mass on s≠0 (WFABC) | Partly — Mathieson et al. test group/period contrasts; Akbari fits one shared time trend with structure covariates only | No — direction estimated once per locus; no concordance statistic across multiple histories | Directly: replace labels with inferred trajectories; the pooled slope becomes per-trajectory slopes plus an across-trajectory replication statistic (our Models 0–3) |
| **2. Admixture-aware selection inference** (adaptive-introgression scans, e.g. Racimo et al. 2015) | Variant (or tract) in an admixed target population, referenced to donor/source populations | Explicit admixture event(s) between named sources; donor sharing via D/f4/S*-type statistics | Background neutral introgression rate | Enrichment of functional/phenotypic variants among introgressed tracts; unusually high-frequency introgressed haplotypes | Yes — per admixed population and per source | No — one admixture history at a time; no across-history replication test | Supplies the logic of conditioning selection inference on a specified mixing event; our trajectories generalize a discrete admixture event to continuous turnover |
| **3. Population-graph methods** (TreeMix, Pickrell & Pritchard 2012; qpGraph, Patterson et al. 2012; landscape population graphs, Dyer & Nason 2004) | Allele-frequency vector of a population; graph nodes = populations, edges = drift/admixture weights | The graph *is* the model: fitted admixture/drift topology with edge weights | None — structure is the target, not the nuisance | Not a selection method; graph-fit residuals can flag selection as model misspecification (exploited by GRoSS) | Edges are population-specific; selection not modeled | No | Graph edges are exactly the "distinct histories" we want as replicates; a temporal version of the graph = our trajectory set; drift weights give the machinery for trajectory covariance (our Model 2) |
| **4. PolyGraph** (Racimo, Berg & Pickrell 2018, Genetics 208:1565) | Set of GWAS trait loci; allele-frequency change along one branch of an admixture graph | Fixed admixture graph fitted from neutral loci | Neutral drift along branches (characterized from the rest of the genome); graph uncertainty | Polygenic allele-frequency change along a branch exceeds the neutral expectation (Q_B statistic) | Yes — branch-specific changes; branches are different populations' histories | Partly — tests one branch at a time; could in principle ask whether a trait rises along several branches, but no such concordance test is implemented | Replace GWAS loci with a single candidate locus and the branch with an inferred trajectory; our concordance statistic is the "same direction along several branches" test PolyGraph never performs |
| **5. GRoSS** (Refoyo-Martínez et al. 2019, Genome Res. 29:1506) | Single locus, tested against an admixture-graph null | Admixture graph; drift along edges characterized from neutral loci | Population history (fully modeled, not removed) | Locus allele frequencies at the tips are improbable under the fitted graph (S_B statistic, branch-specific) | Yes — selection localized to specific branches/populations | No — no formal test of concordant selection on the same locus along multiple branches | Closest single relative: a structure-aware single-locus null. We extend it from cross-sectional (modern tips + static graph) to temporal (dated samples + inferred trajectories) and add explicit replication |
| **6. CLUES / CLUES2** (Stern et al. 2019, PLoS Genet. 15:e1008384; Vaughn & Nielsen 2024, MBE 41:msae156) | Single allele in one (approximately panmictic) population; full allele-frequency trajectory as latent variable | Essentially none — single-population Wright–Fisher diffusion; ancient genotypes and ARGs enter as HMM emissions | Nothing explicit; population structure is *not* modeled (exchangeability assumed) | Full-likelihood LR on selection coefficient s; reconstructed allele history inconsistent with neutrality; CLUES2 allows time-varying s | Not natively — one population at a time | No — per-population s estimates could be compared informally, but no concordance statistic or covariance model | Complementary gold standard: CLUES2 gives locus-specific s and allele history for validating our top loci (AMY, LCT); our statistic is genome-wide and cheap; CLUES can confirm, not scan |
| **7. Ancestry-stratified CLUES** (emerging practice; no distinct published method named "CLES" found — see note above) | Locus within a local-ancestry background of admixed individuals | Local ancestry conditioning: run CLUES within ancestry-specific haplotype panels | Demography within each ancestry background | Selection acting on an ancestry-specific allele copy (s per ancestry background) | Yes — by ancestry background | Yes, in principle — across ancestry backgrounds within one population | Closest conceptual cousin: tests selection *within* distinct backgrounds. Our trajectories generalize "ancestry background" to whole inferred population histories, and add a formal cross-trajectory test |
| **8. Local-ancestry selection** (local ancestry deviation, Bhatia et al. 2014, AJHG; Tractor, Atkinson et al. 2021, Nat. Genet.; aDNA benchmark, Mol. Biol. Evol. 43:msag236, 2026) | Local-ancestry proportion at a chromosomal position within admixed individuals | Mosaic of ancestry tracts from reference populations; genome-wide admixture proportion as baseline | Neutral admixture proportions; drift since admixture; reference-panel error | Deviation of local ancestry from genome-wide expectation (excess ancestry from the population that carried the beneficial allele); Tractor: ancestry-specific association effects | Yes — per admixed population and ancestry source | No — tests selection after one admixture event within one population | Provides within-trajectory direction evidence robust to pooling: does the allele rise *within* each ancestry component? Complements our approach when admixture is recent; ours works when tracts have decayed |
| **9. Simon & Coop drift/migration/selection decomposition** (Simon & Coop 2024, PNAS, doi:10.1073/pnas.2312377121) | Genome-wide variance in allele-frequency change Δp between two dated time points | Explicit drift/gene-flow model between time points; covariances across loci and intervals | Drift and gene flow — they are the quantities being quantified; for locus-level testing they form the null | Genome-wide selection contribution via covariances between time intervals; locus outliers vs the decomposition | Yes — per time-series (applied to two European series) | No — each series analyzed separately; no cross-series concordance statistic | Their decomposition is the right null model for our per-trajectory slopes; their covariance machinery motivates Model 2, and their finding that gene flow dominates Δp motivates trajectory-level analysis |
| **10. Bayenv / LFMM** (Coop et al. 2010; Günther & Coop 2013; Frichot et al. 2013; Caye et al. 2019) | SNP allele frequency across present-day populations vs an environmental variable | Population covariance matrix from control loci (Bayenv) or latent factors (LFMM) | Population structure covariance — the entire point | Allele–environment association beyond structure (Bayenv Bayes factors; LFMM z-scores) | Implicitly — a locus tracks an environment, but no per-population selection parameter | No — no time axis at all | Methodological template: estimate covariance from genome-wide neutral loci and test candidates against it. We transplant this to the time axis (time as the "environment", drift as covariance) |
| **11. Convergent / parallel adaptation** (Lee & Coop 2017, Genetics 207:1591 — DMC; reviewed in Phil. Trans. R. Soc. B 374:20180236) | A shared adaptive allele/phenotype across multiple populations; modes: independent mutation vs shared standing variation vs gene flow | Explicit demographic relationship among the sharing populations (required to distinguish modes) | Shared ancestry and gene flow (they confound "independent" convergence) | Coordinated allele-frequency shifts in the same direction, more than expected given the populations' relatedness | Yes — populations can adopt different modes/alleles | **Yes — core capability**: "surprise" at convergent shifts given covariance; but cross-sectional, not temporal | Their surprise statistic is the static cousin of our concordance statistic; we add the time axis and inferred (unlabeled) histories |
| **12. Phylogenetic comparative methods** (independent contrasts, Felsenstein 1985; PGLS, Grafen 1989; phylogenetic ANOVA, Garland et al. 1993; OU-regime models, Butler & King 2004; SURFACE, Ingram & Mahler 2013; convergence metrics, Stayton 2015) | Species/population trait value at tree tips | Phylogenetic covariance tree (Brownian-motion-derived V matrix) | Shared phylogenetic history (non-independence of tips) | Trait–predictor association robust to phylogenetic covariance; convergence = repeated movement toward the same optimum (regime shifts, C-metrics) | Yes — regime-dependent optima; clade-specific parameters | **Yes** — convergence across independent clades is exactly "repeated directionality across histories"; but requires a known tree and tip values; no allele-frequency time-series semantics | The most direct statistical template for our Models 2–3: trajectory relatedness plays the role of the phylogeny, per-trajectory slope estimates are the tip values, and the effective number of independent trajectories is the effective sample size after V-correction |
| **13. Trajectory inference** (single-cell pseudotime: Monocle, Trapnell et al. 2014; diffusion pseudotime, Haghverdi et al. 2016; RNA velocity, Bergen et al. 2020; aDNA practice: PCA/ADMIXTURE-through-time transects) | Individual cell / individual ancient genome placed on a latent temporal manifold | Low-dimensional manifold or kNN graph connecting individuals; time implicit (pseudotime) or observed | Technical variation (batch, cell cycle ≈ damage, coverage); individual ancestry variation for aDNA | Not a selection method — evidence concerns ordering/lineage commitment (expression change along pseudotime) | No selection semantics | No — trajectories describe history, not locus behavior | Supplies the inference machinery (diffusion on kNN graph, time-aware clustering) our traj-prototype/traj-alt tasks borrow; crucially we must add genetic-descent constraints: a pseudotime ordering is not a demographic trajectory |
| **14. Continuity/replacement detection without archaeological labels** (qpWave/qpAdm, Patterson et al. 2012; qpAdm assessment, Harney et al. 2021; f4 symmetry tests; Legofit, Rogers) | A population (set of dated individuals) vs candidate ancestors/sources | Rank of the allele-frequency covariance of the target against outgroup/reference sets; explicit source sets | Measurement noise | Continuity = later population consistent with drift from earlier one (rank not increased, qpAdm fits); replacement = no fit without additional ancestry streams | Yes — per region and time series | No — purely demographic; no locus-level semantics | qpWave rank tests give a label-free statistical definition of "same trajectory vs new trajectory" — the exact primitive our trajectory segmentation needs (Decision 1); qpAdm weights provide mixture coefficients along a trajectory |

## Key observations

### Structure representation spectrum
Methods represent structure along a spectrum from:

- **Explicit graphs**: TreeMix, qpGraph, PolyGraph, GRoSS (a fully specified demographic model)
- **Covariance matrices**: Bayenv, the Akbari GLMM (statistical approximation)
- **Latent factors**: LFMM (data-driven structure)
- **Conditioning**: ancestry-stratified and local-ancestry methods (ancestry as a conditioning variable)
- **Ignored**: CLUES (panmictic assumption)

### Nuisance treatment
Nearly all selection-scanning methods treat structure as nuisance to be corrected or
conditioned away. Only convergent-adaptation methods (row 11) and ancestry-stratified
practice (row 7) partially use structure as a *replication* framework. This is the
central design tension our project addresses.

### Directionality detection
Only convergent-adaptation methods (row 11) and phylogenetic convergence tests (row 12)
explicitly test for repeated directionality, and both require the independent lineages
(or tree) to be known in advance and act on cross-sectional data. No surveyed method
tests repeated directionality *through time* across *inferred* (unlabeled) population
histories.

## Explicit novelty statement

### What existing methods already solve

1. **Pooled, structure-corrected time-trend scans** — Akbari et al. (2026) GLMM detects
   consistent temporal trends genome-wide while correcting structure as nuisance.
2. **Single-locus, single-population full-likelihood selection inference including aDNA** —
   CLUES/CLUES2, with ARG-based emissions and time-varying s.
3. **Structure-aware selection nulls with an explicit (present-day) admixture graph** —
   GRoSS (single locus, S_B) and PolyGraph (polygenic trait, Q_B).
4. **Genome-wide decomposition of allele-frequency change into gene flow / drift /
   selection** between dated time points — Simon & Coop (2024).
5. **Post-admixture selection via local-ancestry deviation** — Bhatia et al. (2014);
   ancestry-specific association — Tractor (Atkinson et al. 2021).
6. **Convergence-mode inference given known demography, cross-sectionally** — Lee & Coop
   (2017) DMC.
7. **Phylogenetically corrected association and convergence machinery** — IC/PGLS/OU-regime
   methods.
8. **Label-free continuity/replacement primitives** — qpWave/qpAdm rank tests.

### The gap our proposed method fills (what would be newly solved)

No surveyed method jointly:

- (a) infers a small set of temporal population trajectories directly from genome-wide
  aDNA, without archaeological labels;
- (b) estimates per-trajectory directional change at each locus (our Model 0);
- (c) tests *replicated same-direction change across trajectories* with an explicit
  correction for trajectory non-independence (Models 1–2; effective number of
  independent trajectories);
- (d) directly quantifies what a pooled structure correction (the Akbari GLMM)
  suppresses relative to that replication evidence.

Existing methods either treat structure purely as nuisance (rows 1, 10), require the
history to be specified in advance (rows 2–5, 11, 14), are single-population at a time
(rows 6, 7), are cross-sectional rather than temporal (rows 5, 10, 11), or are
demographic-only with no selection semantics (rows 3, 14). The proposed combination —
trajectory-specific slopes, cross-trajectory meta-analysis with covariance correction,
and a concordance statistic with an effective-trajectory-number denominator — is novel
as a *design*, even though each ingredient has a static or single-population precedent.

### What remains open (and is therefore our research risk, not our novelty)

- Whether stable ancestry trajectories can be inferred at all from sparse aDNA
  (**Decision 1**; rows 13–14 supply primitives but no turnkey solution).
- Calibration of the replication statistic under neutral demographic simulations
  (**Decision 2**).
- Power and false-positive trade-offs versus the Akbari statistic (**Decision 3**).
- Quantifying trajectory non-independence before treating trajectories as replicates
  (row 12 supplies the statistical machinery; the population-genetic instantiation is new).
- Verification of what "CLES" refers to in the Akbari supplement (open item).
- Generality beyond West Eurasian aDNA.

## References

- Akbari, A. et al. (2026) Ancient DNA reveals pervasive directional selection across
  West Eurasia. *Nature*, doi:10.1038/s41586-026-10358-1. Preprint: PMC11419161.
- Atkinson, E.G. et al. (2021) Tractor uses local ancestry to enable the inclusion of
  admixed individuals in GWAS and to boost power. *Nature Genetics* 53, 195–204.
- Bergen, V., Lange, M., Peidli, S., Wolf, F.A. & Theis, F.J. (2020) Generalizing RNA
  velocity to transient cell states through dynamical modeling. *Nature Biotechnology*
  (RNA velocity).
- Bhatia, G. et al. (2014) Genome-wide scan of 29,141 African Americans finds no
  evidence of directional selection since admixture. *Am. J. Hum. Genet.*
- Butler, M.A. & King, A.A. (2004) Phylogenetic comparative analysis: a modeling
  approach for adaptive evolution. *Am. Nat.* 164, 683–695 (unverified: pages).
- Caye, K. et al. (2019) LFMM 2: fast and accurate inference of gene–environment
  associations. (LFMM2.)
- Coop, G. et al. (2010) Using environmental correlations to identify loci underlying
  local adaptation. *Genetics* (bayenv).
- Dyer, R.J. & Nason, J.D. (2004) Population graphs: the graph theoretic shape of
  genetic structure. *Molecular Ecology* 13, 1713–1727 (unverified: pages).
- Felsenstein, J. (1985) Phylogenies and the comparative method. *Am. Nat.* 125, 1–15.
- Foll, M., Shim, H. & Jensen, J.D. (2015) WFABC: a Wright–Fisher ABC-based approach for
  inferring effective population sizes and selection coefficients from time-sampled
  data. *Molecular Ecology Resources* 15, 1174–1183 (unverified: pages).
- Garland, T., Harvey, P.H. & Ives, A.R. (1993) Phylogenetic analysis of covariance by
  computer simulation. *Systematic Biology* 42, 265–292 (unverified: pages).
- Grafen, A. (1989) The phylogenetic regression. *Phil. Trans. R. Soc. B*.
- Günther, T. & Coop, G. (2013) Robust identification of local adaptation from allele
  frequencies (bayenv2). *Genetics*.
- Haghverdi, L., Büttner, M., Wolf, F.A., Buettner, F. & Theis, F.J. (2016) Diffusion
  pseudotime robustly reconstructs lineage progression. *Nature Methods* 13, 845–848
  (unverified: pages).
- Harney, E., Patterson, N., Reich, D. & Wakefield, J. (2021) Assessing the performance
  of qpAdm. *Genetics* 217, iyaa045.
- Ingram, T. & Mahler, D.L. (2013) SURFACE: detecting convergent evolution using
  stochastic process models. *Methods in Ecology and Evolution* (unverified: pages).
- Lee, K.M. & Coop, G. (2017) Distinguishing among modes of convergent adaptation using
  population genomic data. *Genetics* 207, 1591–1619.
- Local ancestry inference identifies robust evidence of selection in Neolithic Europe.
  *Mol. Biol. Evol.* 43, msag236 (2026).
- Mathieson, I. et al. (2015) Genome-wide patterns of selection in 230 ancient Eurasians.
  *Nature* 528, 499–503 (unverified: pages).
- Patterson, N. et al. (2012) Ancient admixture in human history. *Genetics* 192,
  1065–1093 (unverified: pages).
- Pickrell, J.K. & Pritchard, J.K. (2012) Inference of population splits and mixtures
  from genome-wide allele frequency data (TreeMix). *PLoS Genetics*.
- Racimo, F., Sankararaman, S., Nielsen, R. & Huerta-Sánchez, E. (2015) Evidence for
  archaic adaptive introgression and its functional importance in humans.
- Racimo, F., Berg, J.J. & Pickrell, J.K. (2018) Detecting polygenic adaptation in
  admixture graphs (PolyGraph). *Genetics* 208, 1565–1584 (unverified: pages).
- Refoyo-Martínez, A. et al. (2019) Identifying loci under positive selection in complex
  population histories (GRoSS). *Genome Research* 29, 1506–1521 (unverified: pages).
- Rogers, A. Legofit: software for estimating admixture proportions and population
  history from allele frequencies (unverified: year/venue).
- Simon, A. & Coop, G. (2024) The contribution of gene flow, selection, and genetic drift
  to five thousand years of human allele frequency change. *PNAS*,
  doi:10.1073/pnas.2312377121.
- Stayton, C.T. (2015) The definition, recognition, and interpretation of convergent
  evolution, and two new measures (C1–C4 metrics). *Evolution* (unverified: pages).
- Stern, A. et al. (2019) An approximate full-likelihood method for inferring selection
  and allele frequency trajectories from DNA sequence data (CLUES). *PLoS Genetics* 15,
  e1008384.
- Trapnell, C. et al. (2014) The dynamics and regulators of cell fate decisions are
  revealed by pseudotemporal ordering of single-cell data (Monocle). *Nature
  Biotechnology*.
- Vaughn, A. & Nielsen, R. (2024) Fast and accurate estimation of selection coefficients
  and allele histories from ancient and modern DNA (CLUES2). *Mol. Biol. Evol.* 41,
  msae156.
