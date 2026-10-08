# Directional consistency: bounded process calibration and AMY feasibility

## Status and preregistration

[observation] Run on 2026-10-08, fixed before results in docs/deme-rebuild-prereg.md. Four timestamped amendments distinguish the computational subpanel, boundary safeguard, later uncertainty/ordered-deme feasibility diagnostics, and a secondary weakest-replicate directional statistic. Outputs: analysis/deme_rebuild/; commands: `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python analysis/deme_rebuild.py` and the same prefix with `analysis/directional_consistency.py`.

[inference] **Partial answer, not a selection discovery.** The directional argument is valid under a specified independence/drift/migration model, not under unrestricted demography. SNP held-out calibration passes the prespecified preliminary gate. AMY replication cannot be evaluated: only ONE sufficiently sampled long CN series, no qualifying sequence of replacing CN-bearing demes, and no independently calibrated CN process. A probability for the AMY ordering is **not identifiable**, not 'nonsignificant.' No biological SNP candidate was interpreted.

## Calibration gate — before locus interpretation

[observation] The fixed, outcome-blind SNP process subpanel comprises demes **25,11,16,15,33,24** (each node >=20 individuals), selected by minimum-node size then span. Of 1,800 disjoint-background SNPs, **212 training and 208 held-out** loci satisfy complete >=8-call/node support and pooled AF [.1,.9]. The support/time design is identical for all controls; frequency scaling supplies coarse AF matching. These are background SNPs, NOT certified neutral loci, and no candidate-specific AF/missingness matching is claimed.

[observation] **Held-out control-locus pass fraction = 13/208 = 0.0625** at two-sided process P<=.05. Exact binomial 95% interval **[0.03370,0.10451]** contains .05. Held-out joint scores have mean **-0.0873**, SD **1.0516**. The preregistered calibration gate **passes**. [inference] This is preliminary consistency with nominal calibration, not proof that the null is correct: LD, background selection and only 208 held-out controls weaken the nominal binomial interval's independence assumption. **CN null: NOT CALIBRATED**; SNP calibration cannot be transferred to CN.

[observation] The fitted standardized increment variance is `a*delta_generations + b*PC_distance_squared + c`, after subtracting conservative pseudo-haploid sampling variance. Fitted **a=0**, **b=0.00101673**, **c=0.000993608**; scaled design condition number **6.345**. [inference] This places **Ne=infinity at the fitting boundary**, NOT evidence that historical drift was zero. Ne is not finitely resolved by these controls; drift and jumps trade off. The centered Gaussian jump approximation is not a donor-resolved admixture model or exact bounded Wright–Fisher simulation. It allows arbitrary signed jumps in variance but does not rule out a directionally biased common donor.

[observation] The separately preregistered training-locus bootstrap (500 joint locus resamples) hits a=0 in **73.8%** of resamples. Drift a quantiles .025/.5/.975 are **0,0,0.000185099**; reciprocal Ne quantiles are **2,701,infinity,infinity**. [inference] These are descriptive fitting uncertainties, not LD-robust demographic confidence limits. A finite independent Ne estimate or explicit donor model is needed before a stronger process claim.

## Independent-unit count and achievable power — before AMY interpretation

[observation] The rebuild retains **42 candidate anchored paths** in **33 fixed geographic cells**, **3,713 assigned people**. They are geographically bounded continuity proxies, not demonstrated genealogical demes. All have **3–5 time nodes**, **35–377 individuals**, **2,151–4,850 years** observed span. No continent-wide blob is retained. Total dated/geolocated/assessment-eligible samples: **18,440**; >=20% called on the independent embedding panel: **12,984**; after biological-ID deduplication: **11,986**. There are **435 retained local nodes** and **269 paths** before duration/size qualification. Only about 31% of coverage/dedup-qualified people join long qualifying paths.

[observation] Standardized background ancestry-profile overlap over **5,746 shared sites** yields maximum pair coefficient **0.7284**, median off-diagonal **0**, and **m^2/sum(R)=15.899** for all 42. On the six SNP-calibration paths, the same ancestry diagnostic is **4.668**. These values discount substantial shared ancestry; they do not establish 16 or 5 independent directional innovations. The inherited estimator clips negative ancestry coefficients and depends on which populations set its reference.

[observation] Training-control temporal slope correlation, shrunk halfway to identity, gives mean-information equivalent **D=6.275** for the six-path test. [inference] This exceeds six because a few estimated correlations are negative; it is NOT a literal number of independent demes or an exponent for a sign probability. Ancestry correlation and correlation of new temporal innovations are different quantities. Actual genealogical independence remains unresolved.

[observation] Conditional Gaussian power for a common per-deme standardized slope effect .25,.5,1,2 is respectively **9.6%,24.0%,70.7%,99.9%** at two-sided .05. [inference] These figures depend on correct fitted variance/correlation and presume homogeneous effects; they are not realized power for ancient CN or detection relative to Akbari. With six truly independent fair signs, perfect two-sided agreement has P=.03125; with only four or five, P=.125 or .0625. The actual CN support is D=1 before any independence discount, so replicated-sign testing has no power here.

## Secondary direct directional test — before AMY interpretation

[observation] Amendment 4 was fixed after the original sum-score results but BEFORE computing this new statistic: `T=max(min(Z),min(-Z))`, the weakest standardized deme slope in either direction. Unlike the original sum, a single large slope cannot drive T. Under the fitted Gaussian process correlation, 100,000 seeded null draws give a 95th percentile **-0.1191**. Held-out control-locus pass fraction **9/208=4.327%**, interval **[1.997%,8.055%]**: the secondary directional calibration gate also **passes provisionally**. Conditional common-effect power .25/.5/1/2 is **7.53/15.27/42.77/90.22%**, lower than the sum test as expected. Monte Carlo power uncertainty is at most about .16 percentage points (one SE); this is still conditional power, not demonstrated biological sensitivity.

[inference] This is a separately pre-run registered but **secondary/exploratory** diagnostic, not a retroactive primary-test replacement. Positive T requires every observed slope to agree; the nominal .05 threshold is negative because iid six-way perfect agreement already has probability only .03125. Thus .05 rejection can include a weak near-zero dissenting slope: report sign counts and require T>0 before describing exact observed unanimity. The held-out fraction checks the model's distribution rather than enforcing five percent by ranking those same controls. No biological SNP candidate is interpreted here.

## AMY1 copy-number observations (not SNP proxies)

[observation] analysis/amy.py's CN loader supplies **537 Ancient/Archaic rows**. Unambiguous AADR token matching resolves **522**, all with actual numeric columns 10/17/18 (BP/latitude/longitude); **424** survive rebuild coverage/geography/date/dedup criteria. **87 carriers** belong to **13 qualifying anchored paths**. No Akbari linkage is needed for this analysis; the previous zero-ENA-overlap conclusion was not a barrier to AADR-based trajectory labeling.

[observation] Only **deme 41** has >=3 nodes with >=5 CN carriers and an AMY-only span >=2,000 years: **29 carriers in four nodes**, with supported-node counts **5,10,10** (one additional late node has 4). AMY observed span **3,587 years**; supported-carrier span **3,487 years**. This path is geographically in **55–60N,10–15E**, has **129 total AADR people**, five AADR nodes, and total date span **575–4,686 BP**.

[observation] The supported CN node means (old to recent) are **7.6** at mean **4,258.6 BP** (n=5), **5.9** at **3,714 BP** (n=10), **6.1** at **1,047.5 BP** (n=10). Supported-node slope = **-0.1883 copies/1,000 forward years**. Including the fourth small recent node (mean CN=7.5, n=4, BP=900) gives a descriptive slope **-0.0292**. [inference] This is not a replicated increase, but one noisy, unequally sampled regional series cannot falsify the author's functional-selection claim. The chronology skips much of the starch-agriculture transition and is not an independent Eurasia/Andes comparison.

[observation] For the BETWEEN-successive-demes reading, **five paths** have >=5 AMY carriers; no cell contains more than **two** such paths. **Zero qualifying nonoverlapping successions of >=3 CN-bearing demes** meet the additional preregistered feasibility definition. [inference] There is no honest ordered-deme tail probability to estimate here. Treating unrelated geographic paths' mean dates as a succession would reinstate the defect being repaired.

[inference] **AMY drift-plus-jump ordering probability: NOT IDENTIFIABLE.** CN mutation/haplotype copy-count variance, donor distribution and measurement error are not calibrated by neutral SNPs, and only one long CN series is supported. Do not report 2^(-D), a pooled CN slope, or a SNP-calibrated P as evidence for the accepted AMY argument.

## Validation, deviations, and limits

[observation] Eleven source-named regressions pass in analysis/test_deme_rebuild.py: exact high-first encoding and arbitrary site order; actual cross-bundle call concordance; adjacent-bin duration; geographic separation; gap splitting; anchor constraint; 8,000-year cap; full cumulative random-walk covariance; forward-time direction/shared-innovation discount; a neutral common-donor counterexample with six concordant slopes; and a weakest-replicate score that cannot be driven by one outlier. The first rebuild run failed on integer division output in the inherited relatedness helper; the caller now supplies floats. The first process run failed on boundary control means, leading to preregistered Amendment 2 before any P-values. No control distribution was remapped to manufacture a 5% rate.

[observation] The inherited load_aadr.py decoder used low-first bits; direct high-first decoding agrees at **339,933/339,933** paired called HO/1240K values, versus **178,989/260,356=68.75%** low-first. The corrected loader is covered by real-data and byte-level regression tests. [inference] Earlier analyses using this loader have site identity/missingness/region-exclusion errors; those numerical results require rerunning, especially relatedness. trajectory_v2.py already uses high-first encoding, so its entire output is not automatically invalidated.

[inference] Remaining limitations: hand-specified geographic/time/genetic gates; no continuity bootstrap or qpAdm donor decomposition; no close-family exclusion beyond biological duplicate removal; low-coverage mean-imputed PCA and sparse nodes; genome-wide LD/ascertainment; Gaussian boundary approximation; unresolved finite Ne; centered rather than directional donor jumps; limited SNP panel; no CN-null calibration or Andes replication. These make the project a feasibility result, not a completed causal selection test.

[observation] No raw data, Lean or existing authors' conclusions were overwritten. No git command was run; parent owns review/commit. Provenance manifests disable git probing and record code/prereg/input metadata/output hashes; large raw genotype files are recorded by path and size, not rehashed during the bounded session.

[inference] Recommended next analysis: prespecified richer independent population series with donor-resolved admixture and a CN haplotype transition model, ideally including the Andes data. First rerun inherited loader-dependent analyses. Do not use this session to claim either AMY selection proved or AMY selection refuted.

## All qualifying rebuilt paths — audit table

[observation] Every row below is measured from the fixed rebuild. Cells show latitude/longitude bounds in degrees, not inferred homeland. Sizes, BP intervals and nodes are total AADR support; they do NOT imply AMY support or independence. Full ancestry matrix is analysis/deme_rebuild/ancestry_overlap.tsv.

| Deme | Latitude | Longitude | n | nodes | min node n | BP min–max | span yr |
|---:|---|---|---:|---:|---:|---|---:|
| 0 | 20–25 | 100–105 | 50 | 3 | 8 | 1479–3808 | 2329 |
| 1 | 30–35 | 35–40 | 64 | 4 | 9 | 698–3775 | 3077 |
| 2 | 30–35 | 110–115 | 72 | 5 | 8 | 2100–6950 | 4850 |
| 3 | 35–40 | -10–-5 | 43 | 3 | 8 | 3546–5950 | 2404 |
| 4 | 35–40 | -5–0 | 53 | 4 | 9 | 400–3981 | 3581 |
| 5 | 35–40 | 10–15 | 60 | 3 | 8 | 2238–4725 | 2487 |
| 6 | 35–40 | 25–30 | 49 | 3 | 10 | 399–2550 | 2151 |
| 7 | 35–40 | 35–40 | 57 | 3 | 13 | 3337–5767 | 2430 |
| 8 | 40–45 | -5–0 | 56 | 3 | 11 | 3275–5700 | 2425 |
| 9 | 40–45 | 0–5 | 83 | 5 | 10 | 1038–5867 | 4829 |
| 10 | 40–45 | 5–10 | 45 | 3 | 9 | 2250–4415 | 2165 |
| 11 | 40–45 | 10–15 | 96 | 3 | 27 | 300–2837 | 2537 |
| 12 | 40–45 | 15–20 | 49 | 3 | 12 | 1029–3631 | 2602 |
| 13 | 40–45 | 20–25 | 59 | 3 | 16 | 246–2900 | 2654 |
| 14 | 40–45 | 25–30 | 58 | 3 | 9 | 4051–6833 | 2782 |
| 15 | 40–45 | 25–30 | 74 | 3 | 22 | 6229–8563 | 2334 |
| 16 | 40–45 | 40–45 | 87 | 3 | 22 | 1050–3945 | 2895 |
| 17 | 40–45 | 40–45 | 93 | 5 | 9 | 2256–6850 | 4594 |
| 18 | 40–45 | 45–50 | 51 | 3 | 11 | 2022–4850 | 2828 |
| 19 | 45–50 | 0–5 | 35 | 3 | 10 | 4030–6785 | 2755 |
| 20 | 45–50 | 5–10 | 79 | 4 | 15 | 1350–4635 | 3285 |
| 21 | 45–50 | 10–15 | 65 | 3 | 12 | 1180–3962 | 2782 |
| 22 | 45–50 | 15–20 | 377 | 5 | 17 | 702–4718 | 4016 |
| 23 | 45–50 | 15–20 | 86 | 3 | 12 | 5019–7725 | 2706 |
| 24 | 45–50 | 20–25 | 279 | 4 | 21 | 550–3966 | 3416 |
| 25 | 45–50 | 20–25 | 124 | 3 | 30 | 5097–7340 | 2243 |
| 26 | 45–50 | 20–25 | 73 | 3 | 20 | 5100–7700 | 2600 |
| 27 | 45–50 | 25–30 | 89 | 4 | 12 | 2194–5967 | 3773 |
| 28 | 45–50 | 30–35 | 45 | 4 | 9 | 4103–7943 | 3840 |
| 29 | 45–50 | 100–105 | 71 | 4 | 9 | 540–3632 | 3092 |
| 30 | 50–55 | -5–0 | 174 | 5 | 12 | 856–4939 | 4083 |
| 31 | 50–55 | 0–5 | 113 | 4 | 19 | 129–3994 | 3865 |
| 32 | 50–55 | 0–5 | 123 | 3 | 17 | 186–2912 | 2726 |
| 33 | 50–55 | 10–15 | 163 | 5 | 21 | 199–4768 | 4569 |
| 34 | 50–55 | 10–15 | 62 | 3 | 13 | 1000–3945 | 2945 |
| 35 | 50–55 | 10–15 | 111 | 3 | 21 | 2197–4719 | 2522 |
| 36 | 50–55 | 10–15 | 107 | 4 | 18 | 4200–7200 | 3000 |
| 37 | 50–55 | 50–55 | 59 | 3 | 12 | 3475–5964 | 2489 |
| 38 | 50–55 | 55–60 | 38 | 3 | 8 | 2051–4850 | 2799 |
| 39 | 50–55 | 90–95 | 66 | 3 | 11 | 2219–4937 | 2718 |
| 40 | 55–60 | -5–0 | 46 | 3 | 13 | 1000–3993 | 2993 |
| 41 | 55–60 | 10–15 | 129 | 5 | 11 | 575–4686 | 4111 |

[observation] Additional implementation deviation from the embedding wording: the actual embedding panel excludes ALL of chromosomes 1 and 6, not only AMY/MHC intervals (the held-out control panel was prespecified to exclude both whole chromosomes). [inference] This conservative exclusion avoids target-region leakage but narrows the genome background; it is disclosed rather than presented as exact adherence to the region-only embedding wording. No panel was reselected after results.
