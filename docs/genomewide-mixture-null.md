# Genome-wide measured-mixture residual validation

[observation] Second-pass session 2026-10-08. Protocol and pre-run amendments: docs/mixture-null-prereg.md. Implementation: analysis/genomewide_mixture_null.py. All scans below are **genome-wide-distributed control panels**, NOT full 1.24M/52M-variant selection discovery. No biological SNP was interpreted. All weights are measured-label one-hot indicators, not fractional ancestry; within-label turnover remains unmodeled.

## Detection floors BEFORE locus interpretation

[observation] With qualifying public AADR ancestry/period labels, the approximate common-ramp 80% floor is **.1134 allele-frequency/kyr** (1240K median; 10–90% .0805–.1431) and **.1131/kyr** (HO; .0830–.1450). The corresponding pgen panel has **zero supported loci**: its floor is not identifiable, not zero. It has only 115 exact shared individuals before locus support restrictions.

[observation] With independently checked longitudinal-deme labels and the TRAIN-only variance correction specified below, median 80% floors on HELD OUT are **.0309 frequency/kyr AADR** (10–90% .0233–.0424) and **.0260/kyr pgen** (.0197–.0291). For example, these are approximately 12.4 and 10.4 percentage points over 4kyr, respectively. Raw floors .0286/.0210 were underestimated relative to temporal background variation and MUST NOT be used to interpret loci because their raw null fails.

[inference] These floors test a common signed within-group linear frequency ramp; they do not guarantee separate significance in each deme or a selection coefficient s. Large enough ancestry-adjusted temporal excesses CAN pass this bar; short series and weak changes cannot. Broad synchronized shifts that also characterize TRAIN controls are partly included in the measured null, limiting power to distinguish very widespread selection from other genome-wide processes.

## Measured inputs and orientation

[observation] AADR: data/raw/aadr-v66.p1/v66.p1_{1240K,HO}.aadr.patch.PUB.{geno,ind,snp}; each bundle's corresponding v66.p1_{1240K,HO}.aadr.PUB.anno supplies BP dates, biological IDs and quality assessments. Genotypes from analysis/load_aadr.py are **ALT dosage**, missing code3 is excluded, not counted as an allele. Divide ALT dosage by2 for individual frequency observations, not two independent pseudo-haploid chromosomes. No binomial diploid sample-size fiction is used.

[observation] Orientation independently rerun on 2000 exact-allele matched SNPs and 300 shared individuals: low-first corr **-.1168**, high-first reference dosage **-.8103**, high-first inverted ALT **+.8103**, agreement **.7223**. analysis/test_aadr_orientation.py passes. The existing encoding/deme tests also pass **11/11**, including actual matching HO/1240K calls. Their first invocation without PYTHONPATH failed import; rerun with PYTHONPATH=. passed.

[observation] Akbari pgen: data/raw/akbari/tier2/AncientDNA_plus_EUR1KG_AllVariants_NoQC_GRCh37_autosomes.{pgen,pvar,psam}; analysis/load_akbari_geno.py returns allele_idx=1 ALT dosage. This scan additionally excludes nonfinite/negative/>2 values, because sentinel values cannot be trusted as allele observations. Samples match exact PSAM IDs, variants match CHROM/POS/REF/ALT exactly; no strand or coordinate guessing.

## Ancestry/period-label feasibility

[observation] Ancient 0<BP<=14000, no FAIL/CRITICAL, biological ID deduplication gives **3512** observed labels in each AADR bundle. Only **17** satisfy n>=16, >=3 distinct dates and >=2000yr span: **506** people. Public labels include China_MLBA, Cuba_PreCeramic, Czechia_N, France_Mesolithic, France_Occitanie_Gard_AvendelaBoucle_N, Italy_Sardinia_N, Kenya_PastoralN, Russia_AltaiKrai_N, Russia_Irkutsk_EBA, Russia_MBA_Caucasus, Russia_Vologda_Mesolithic, Serbia_IronGates_Mesolithic, Switzerland_EBA-2, Turkey_C, Turkey_N, Turkey_PPN, Ukraine_N. This is a worldwide feasibility panel, not a western-Europe-only population history.

[inference] Period-specific labels mostly fail temporal support. AADR can nevertheless test some broad-label changes, but genotype intersection destroys usable pgen support in this panel. A negative result from this label choice would NOT test most lineages or AMY. Qualifying a label by time depth does not prove genetic continuity.

## Every control-locus pass fraction

[observation] Controls: every200th eligible autosomal .snp index, chromosomes1/6 excluded; choose evenly spaced1200 across all exact pvar matches. No slope/X filtering. Locus gate: >=100 supported called people; >=2 groups each >=8 calls, >=3 dates and >=2kyr span; AF .1–.9. TRAIN first half, HELD OUT second half in fixed genomic input order. Raw HC3 uses measured member means, centered within-group time, df=N-K-1. Calibration gate requires held-out N>=100 and exact 95% binomial interval containing .05.

| Panel / null | Eligible controls | Held out | P<=.05 | Fraction | Exact 95% CI | Gate |
|---|---:|---:|---:|---:|---:|---|
| Public labels AADR1240K / HC3 | 481 | 241 | 15 | 6.22% | 3.53–10.06% | provisional pass |
| Public labels AADR HO / HC3 | 476 | 238 | 17 | 7.14% | 4.22–11.19% | provisional pass |
| Public labels pgen / HC3 | 0 | 0 | undefined | undefined | undefined | insufficient support |
| Longitudinal demes AADR / raw HC3 | 790 | 395 | 32 | 8.10% | 5.61–11.24% | **FAIL** |
| Longitudinal demes pgen / raw HC3 | 779 | 390 | 31 | 7.95% | 5.46–11.09% | **FAIL** |
| Longitudinal demes AADR / TRAIN variance | 790 | 395 | 24 | 6.08% | 3.93–8.91% | provisional pass |
| Longitudinal demes pgen / TRAIN variance | 779 | 390 | 15 | 3.85% | 2.17–6.26% | provisional pass |

[observation] These are raw nominal P-values, NOT held-out ranks engineered to return5%. Every tested control panel covers chromosomes2–5 and7–22. AADR/HO and pgen share people/alleles; their agreement is not independent biological replication. Confidence intervals assume independent loci, while LD can reduce effective control N; every200th spacing is not a full LD-pruning guarantee.

[inference] Both raw longitudinal nulls were broken and were NOT used to interpret a locus. The dated repair constrains a single extra variance parameter from different loci. Passing the finite-panel gate does not prove exact5% calibration; especially 7.14% with a wide interval remains imprecise. No discovery claim follows.

## Longitudinal fallback validated, not blindly trusted

[observation] Source analysis/deme_rebuild.py builds 5° geographic-cell, adjacent1kyr-bin, eight-PC ancestry paths; nearest one-to-one matching requires centroid distance <=1.5×median cell radius, anchor distance <=2×radius and total span<=8kyr. Embedding is6000 fixed-seed background SNPs excluding chr1/6; no AMY information is used. Our validation checked every saved row ID, biological identity, numeric date, deduplication, node count/date average, adjacent temporal bins and geographic-cell membership against current AADR. All check. Source constraints are regression-tested; the large original PCA rebuild was NOT rerun in this session, and historical genealogical continuity remains inferred, not observed.

[observation] Independently checked artifacts yield **42** qualifying paths, **3713** people, **149** nodes, **35–377** people/path and **2151–4850yr** actual spans. **2265** people match pgen exactly. **Two** controls overlap the embedding panel and are removed BEFORE residual testing. SNP-specific support varies; AADR median n2324/K41, pgen n2265/K39.

[inference] These have substantially more longitudinal information than public period labels but are candidate paths, not guaranteed independent genealogical innovations. The common slope is an excess-average test, NOT proof of unanimous replicated selection. Shared drift/migration within these groups is precisely why individual HC3 alone fails.

## One additional measured genome-wide parameter, not locus-specific donors

[inference] For bundle b, estimate on TRAIN loci only:

    tau_b² = max(0, mean_train[(beta_l² - SE_l²)/(AF_l(1-AF_l))])
    V_l = SE_l² + tau_b² AF_l(1-AF_l)
    test beta_l / sqrt(V_l), Student-t df=L_train-1.

[inference] Null mean is fixed zero; no donor frequency, locus slope correction or mean centering is fitted. This is a background temporal-slope second moment, not an Ne estimate. TRAIN loci's own K component means are measured and charged to their slope df. Each held-out locus's measured baseline means cost K_l and its common slope costs1; their constant uncertainty cancels in the centered numerator. Genome-wide correction costs **one** TRAIN variance parameter per bundle (395/389 training loci), **zero** fitted correction parameters on validation loci. HC3 residual second moments are sampling-uncertainty estimates, not adjustable demographic knobs.

[observation] AADR tau²=**.00010004**; TRAIN bootstrap95% range **0–.00023486**. Pgen tau²=**.00015338**, range **.00009288–.00022081** (500 seeded resamplings). Variance is not robustly nonzero in the AADR training panel. No new correction was selected after observing repaired held-out results.

[inference] Ancestry labels/paths themselves are independently measured covariates, NOT uncharged locus-fitted demographic parameters. The background inference has many numerical estimates (background site means/scales, eight PCA loadings/scores, k-means centroids/radii and fixed linking rules), so pretending its complexity is zero would be wrong. Its charge is **locus holdout**: no tested site enters the saved6000-site embedding or parameter estimation. Up to6000 means/scales and8×6000 loadings (with orthogonality constraints),8 scores/person, centroid/radius per original node, and cell-level radius summaries are background estimates; they are not identified anew using the test site's allele trajectory. Only149 qualifying nodes are used here; original source reports435 total nodes. The full ancestry estimator is conditioned on, not granted freedom to reshape each candidate trajectory. Weight indicators are observed memberships, not individually fitted admixture weights. This is the held-out-locus form of DOF accounting, not subtracting background PCA dimensionality from each SNP's N.

## Conditional positive-effect power

[observation] Under the registered corrected Gaussian slope approximation, 10000 seeded score-noise replicates across held-out designs give:

| Common frequency ramp/kyr | AADR mean power | Pgen mean power |
|---:|---:|---:|
| 0 | 4.93% | 4.94% |
| .01 | 15.6% | 21.0% |
| .02 | 45.7% | 61.6% |
| .03 | 75.6% | 90.4% |
| .04 | 91.4% | 98.8% |
| .05 | 97.2% | 99.9% |
| .075 | 99.8% | ~100% |
| .1 | ~100% | 100% |

[inference] These conditional score simulations show the variance repair does not annihilate signal, but are NOT independent validation of genotype evolution. The smallest tested-grid effect giving average>=80% power is **.04/kyr AADR**, **.03/kyr pgen**; individual-locus median floors .0309/.0260 differ because averaging heterogeneous power curves differs from power at a median design. CN's separate full within-group dataset injection and pure-composition specificity experiment are reported in docs/amy-mixture-null-results.md.

## Akbari statistics: independent units check, not a neutrality oracle

[observation] analysis/load_akbari.py streams the raw release and exact-joins **1242/1308** unique eligible control IDs across all panels. Median signed X=**.1926755**; **zero** have |X|>5.45; **43** have |X|>=3.16. X was **not centered**. |X|>5.45 is the paper's pi>99% threshold; |X|≈3.16 is approximately FDR50%, not conventional genome-wide significance. Background controls were never selected using these values. Initial public-label-only join had915/956 matches, median X=.190754,32 at>=3.16; after fallback the combined join above supersedes its scope, not its numbers.

[inference] Low X does not prove neutrality, and a background second moment can absorb real broadly distributed selection. Conversely the null is not allowed to invent a high-frequency AMY donor to explain a tested CN trend. A defensible correction predicts composition from independently measured sources and evaluates within-source excess against independent background variation, with explicit power and finite uncertainty.

## Reproducibility and limits

[observation] Commands: .venv/bin/python analysis/genomewide_mixture_null.py; --demes; --calibrate-demes; --join-stats. Outputs, per-control SNP betas/SEs/P-values/support/floors/partitions and JSON summaries are under analysis/mixture_null/. Dated preregistration preceded each added statistic. Seeds are20261008. No git command or raw-data edit occurred.

[inference] This is a bounded validation implementation. No full scan, locus-level FDR, independence-of-trajectories proof, mechanistic selection-coefficient detection floor, fractional ancestry or CN evolutionary calibration was produced. The corrected held-out panels have already been inspected once under the failed raw null; an untouched chromosome/block-pruned external panel is the next validation step, not a retrospective claim of pristine validation.

[observation] Additional new source-named genome-wide regressions pass **4/4**: direct public AADR ALT decoding, per-locus support/negative missing-sentinel exclusion, frozen TRAIN-only variance invariance to arbitrarily changed held-out slopes, and the tau²=0 boundary. These supplement the5 new CN residual regressions and11 existing source/encoding/path/process tests.
