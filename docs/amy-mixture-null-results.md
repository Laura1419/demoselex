# AMY measured-mixture audit — incremental results

[observation] Appended-session analysis 2026-10-08; procedure registered before runs in docs/mixture-null-prereg.md. New source-derived tables: analysis/mixture_null/. Figure: fig/amy_mixture_null.png. This is the categorical, one-hot **COARSE LIMIT**, not a full fractional-ancestry correction or a calibrated evolutionary CN-process null.

## Detection floor FIRST

[observation] The common within-label slope's observed HC3 80%-power approximation is **0.335 copies/kyr**; homoskedastic SE=.155 implies approximately **.433 copies/kyr**. The actual registered heteroskedastic simulation is more conservative: power .4=76.3% (95% CI 73.5–78.9%), .45=85.0% (82.6–87.2%). Thus the refined tested-grid **80% detection floor is .45 copies/kyr**, with crossing between .4 and .45. This corresponds to approximately 1.8 extra copies over 4kyr if that common ramp operates throughout. It is NOT a minimum detectable selection coefficient s.

[observation] For individual labels, approximate observed HC3 80% floors (copies/kyr) are Early farmers .847; WHG .651; EHG .528; Neolithic farmers 2.174; Steppe 3.017; Bronze 2.314; Iron-to-modern 2.312. A singleton CHG has **no identifiable slope**. Floors are conditional quantitative-trait sampling approximations, not power against arbitrary drift/turnover. The low common floor depends on a common signed ramp across components; it is not power for repeated separate discoveries.

[inference] These data CAN test a sufficiently large sustained common change within the coarse measured labels under the specified noise null. They cannot in themselves test s≈.02, prove uninterrupted lineage continuity, or exclude fractional ancestry change inside a label. The label-specific detection floors are often several copies/kyr. Failure to reject is not a refutation of selection.

## Mandatory falsifiability AND specificity

[observation] Fixed seed 20261008; 1000 replicates/effect; actual observed dates, component baselines and measured residual SDs. Pure composition with FLAT within-group CN: **41/1000=4.1% rejection**, exact CI **2.96–5.52%** (contains 5%). The uncorrected pooled trend falsely rejects this same null **99.9%** of the time. No invisible donor was introduced.

| Injected within-label copies/kyr | Residual-test power |
|---:|---:|
| 0 | 4.1% |
| .05 | 6.0% |
| .1 | 9.7% |
| .2 | 25.0% |
| .3 | 51.6% |
| .35 | 64.0% |
| .4 | 76.3% |
| .45 | 85.0% |
| .5 | 91.8% |
| .75 | 99.9% |
| 1, 1.5, 2 | 100% |

[observation] The table reports registered simulated observations, not estimated historical selection. Full exact intervals are in amy_power.tsv and amy_power_refined.tsv.

[inference] The procedure demonstrably rejects true within-label ramps without rejecting pure replacement more than its nominal rate in this conditional experiment. Therefore strong signals are NOT all corrected away. A procedure unable to reject injected positive effects would be worthless for that class; this one passes that necessary but not sufficient validation. Gaussian simulations do not validate biological CN drift or sequencing batches.

## Independent recomputation of supplied findings

[observation] Exact sample join yields **533** Ancient individuals with p2/date/CN. Supplied seven group sizes/means are reproduced, but they sum to **532**, not 533: the release also includes **CHG n=1, mean AMY1=8, age=9724 BP**. There are **eight** labels. Excluding CHG, age-ordered means reproduce **5.45→5.75→5.39→7.15→8.25→6.63→7.25**, with the Bronze drop. Including CHG inserts 8.00 at the oldest position. The panel is dominated by Northern Europe (333); regional labels include 6 NorthAsia and 1 CentralAsia, so “West Eurasian” is an approximate coverage description, not an explicit regional filter.

| Component | n | Mean CN | Unique dates | Span kyr | Slope copies/kyr | HC3 z | P |
|---|---:|---:|---:|---:|---:|---:|---:|
| CHG | 1 | 8.00 | 1 | 0 | undefined | undefined | undefined |
| Early farmers | 22 | 5.455 | 22 | 5.439 | .407 | 1.348 | .193 |
| WHG | 32 | 5.750 | 31 | 5.803 | .043 | .184 | .855 |
| EHG | 49 | 5.388 | 49 | 6.819 | .149 | .789 | .434 |
| Neolithic farmers | 34 | 7.147 | 33 | 3.716 | -.339 | -.437 | .665 |
| Steppe pastoralists | 16 | 8.250 | 15 | 2.921 | .061 | .057 | .955 |
| Bronze age | 46 | 6.630 | 43 | 1.696 | .433 | .524 | .603 |
| Iron to early modern | 333 | 7.249 | 51 | 1.550 | -.057 | -.069 | .945 |

[observation] No supported within-label slope rejects zero. The supplied “EVERY component |z|<=1.3” needs qualification: CHG slope is undefined; registered HC3 gives Early farmers **|z|=1.348**, marginally above 1.3. With ordinary homoskedastic linear-regression SE, its z=**1.279**, reproducing the supplied bound. These are different uncertainty conventions, not a biological contradiction.

[observation] Five coarse old components satisfy >=16 people, >=3 dates and >=2kyr span (Early farmers, WHG, EHG, Neolithic farmers, Steppe). Bronze and Iron do not meet the span threshold despite 46 and 333 observations. The release has dated members over long intervals, NOT repeated sampling of the same individuals, and the labels are ancestry/period categories, not demonstrated persistent demes.

## Degree-of-freedom stress test

[observation] All rows use N=533. K counts actual measured component baselines, with slope df=N-K-1; residual variance is additional uncertainty estimation, not an unobserved donor.

| Fixed partition | K | Sxx (kyr²) | beta copies/kyr | HC3 SE | P | Approx 80% floor |
|---|---:|---:|---:|---:|---:|---:|
| Pooled | 1 | 4031.066 | .218 | .039 | 3.31e-8 | .109 |
| p2 | 8 | 291.971 | .151 | .119 | .206 | .335 |
| p2×2kyr bins | 21 | 39.642 | .609 | .407 | .135 | 1.139 |
| p2×1kyr bins | 34 | 16.065 | .355 | .626 | .571 | 1.753 |
| p2×.5kyr bins | 55 | 5.161 | .168 | 1.166 | .886 | 3.266 |
| p2×.25kyr bins | 81 | 1.823 | -1.381 | 1.967 | .483 | 5.511 |
| Singleton | 533 | 0 | not identified | — | — | infinite/unidentified |

[observation] Seeded half-member measured baselines change the residual beta by **3.89e-16 copies/kyr**. Within-group time centering makes the slope invariant to every constant group baseline, so estimating a component mean cannot cancel an injected within-group slope. Subdividing groups over time is different: it removes the contrasts themselves.

[inference] Exactly at **K=N=533**, unrestricted individual baselines can reproduce every observation and beta is unidentifiable. Unfalsifiability is not merely “more conservative P-values”: the design loses its target information. Already K=81 leaves only .045% of pooled time information and has an approximate floor >5 copies/kyr. These subdivisions are an overflexibility stress test, NOT approved ancestry correction.

## Interpretation (after floors and validation)

[observation] A pooled positive slope (.218 copies/kyr, P=3.31e-8) becomes an imprecise positive within-label slope (.151±.119 HC3 SE, P=.206). [inference] The pattern is **compatible with** composition of measured category means; it does not REQUIRE detectable within-label change. This is not proof that the rise was purely demographic. The observed positive residual, weak power for .1–.2 copies/kyr and unmeasured within-label turnover prevent that stronger conclusion.

## Haplotype source audit

[observation] All **10** Ancient H1/H2 structural labels also occur in Modern rows; Modern has 11 labels (additional H3A3B3). Slot-weighted **diploid carrier** AMY1 means reproduce H1^a **4.299**, H3^r **6.487**, H5 **9.015**, H7 **11.161**, H9 **13.952**. These are NOT haploid structural copy counts; repeated slots weight homozygotes twice.

[observation] Modern H5 slot frequencies: EA **.4360** (1438 slots), WEA **.1864** (1738 slots), PEL **.4471** (170 slots), Pima **.7143** (28 slots). Thus “PEL .46” is not reproduced: the release value rounds to **.45**, not .46. “Europe .18” is approximate: WEA is **.1864**, rounds to .19 and includes Western Asia; it is not an exclusively European bin. EA ~.44 and Pima ~.71 are reproduced. Full population values are in h5_by_p1.tsv, macroregions in h5_by_p2.tsv.

[inference] Shared structural labels show that these data do not establish an independently arisen AMY mutation on each continent. They do NOT prove that all copies with a label are identical by descent, nor that Andean H5 came from a particular East Asian ancestor. Shared labels alone cannot quantify how much Andean elevation is explained or rule out independent selection on standing inherited variation. That stronger supplied causal statement is **not verified** by this release-summary test.

[observation] Further registered descriptive decomposition and longitudinal-deme feasibility will be appended as computed; no pending result is used above.

## Additional registered outputs — measured decomposition and strict series

[observation] **Strict-series detection floor FIRST:** Across all87 carriers in13 continuity-constrained paths, a descriptive common-slope HC3 floor is **1.003 copies/kyr**, much worse than the coarse p2 floor. The sole strict usable path's floor is **1.188 copies/kyr**, n=25 supported carriers over **3.487kyr**. These are uncalibrated CN sampling floors, not process-calibrated significance. They precede the strict-series interpretation below.

[observation] Exact-token independent linkage of533 release individuals to saved, independently checked rebuilt AADR identities retains **424**, with **zero ambiguous tokens**. **87** fall in13 qualifying geographic paths. Only **ONE** path satisfies >=3 nodes, >=5 CN carriers/node, >=2kyr span; it has29 total carriers across4 nodes,25 across3 supported nodes. Thus **zero pairs of adequately supported strict CN time series** are available. Path-level counts/spans are in amy_strict_demes.tsv; supported n and span for the others are mostly5–17 people and .075–.589kyr, or no supported nodes.

[observation] The one strict path's individual-date regression slope is **-.182 copies/kyr**, HC3 SE **.424**, P=.672; all87 carriers' coarse shared within-path slope is **-.058±.358**, P=.871. This differs slightly from the prior lane's node-mean-time slope -.188 because this run uses each individual's AADR date rather than weighted node-average dates. Archive dates are retained for the PRIMARY p2 analysis; AADR dates are used here to respect existing node construction.

[inference] The paper's release has enough time depth/n to test very large changes across broad p2 labels, but **does NOT supply replicated continuity-constrained CN longitudinal series** at the registered support bar. Therefore this dataset cannot test the strong independent-within-lineage AMY replication claim in that strict sense. Its lack of strict replication is a data limitation, not evidence that AMY selection did not occur.

[observation] Measured p2-mean composition prediction has pooled slope **.20725** versus observed **.21822**, a descriptive **94.97%** of the pooled slope. Residual contribution is .01097, because only **7.24%** of pooled time information survives p2 centering. With baselines measured in seeded50%-member subsets, composition slope is **.15024**; prediction RMSE on267 held-out members is **2.759 copies**. Thus the nearly95% in-sample fraction is NOT a precise independently estimated causal contribution; sampling the small components changes it substantially.

[inference] This in-sample decomposition is partly an algebraic consequence of the measured group means and design. It shows compatibility with composition, NOT proof of pure replacement. The exact baseline-invariance of the within-group slope is the inferential advantage: independently measured baseline uncertainty cannot be used to dial away a residual ramp.

[observation] Modern **Karitiana+Surui** H5 frequency is exactly **.2000** across40 haplotype slots, reproducing the supplied Amazonian summary for this specifically named subset. Karitiana alone .2917 (24 slots), Surui .0625 (16); a catch-all “Amazonian population” claim would hide substantial heterogeneity.

## Direct power loss when nuisance baselines proliferate

[observation] The SAME simulated CN matrices are tested with each partition, not regenerated to fit the correction. This controls everything except nuisance-partition flexibility:

| K / partition | Pure composition rejection | Power beta=.5 copies/kyr | Power beta=2 |
|---|---:|---:|---:|
| 1 / pooled | 99.9% | 100% | 100% |
| 8 / p2 | 4.1% | 91.8% | 100% |
| 21 / p2×2kyr | 5.0% | 23.1% | 99.8% |
| 34 / p2×1kyr | 4.4% | 10.4% | 85.4% |
| 55 / p2×.5kyr | 4.2% | 5.4% | 37.1% |
| 81 / p2×.25kyr | 3.0% | 3.8% | 14.1% |
| 533 / singleton | undefined | undefined | undefined |

[inference] At K55 the procedure has almost no power for an effect that K8 detects91.8% of the time. At K81 even a2-copy/kyr ramp is detected only14.1%. These excessively subdivided corrections are practically worthless for the .5 effect; singleton correction is mathematically unfalsifiable. Adding parameters must therefore be justified by independent ancestry measurements and accompanied by loss-of-information/power accounting, not excused solely because the resulting trend is nonsignificant.

## Provenance of the categorical labels

[observation] The extracted source graph_genotyping/assess_graph_genotypes/selection_analysis.Rmd lines160–185 joins master_table_wide_ar.xlsx annotations p2_l0/p2_L1, includes manual category corrections (e.g. NEO816 EarlyFarmers), and maps EarlyFarmers→Early farmers, Europen Farmers→Neolithic farmers, Steppe→Steppe pastoralists, PostNeolithic_BA→Bronze age, EHG&CHG→EHG/CHG. Supplement metadata exports p2 as p2_l0 at lines284–296. This source trace is not a documented measured fractional-ancestry estimator.

[inference] Accordingly, the label-mixture test is the available coarse measurement limit; treating these bins as quantified ancestry fractions or as independently evolving populations would exceed the data. The correct negative statement is “no required detectable common residual under coarse labels,” not “s=.02 was explained away by measured demography.”
