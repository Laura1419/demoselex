# H. AMY analysis

> Owner task: locus-amy. Status: descriptive baseline complete; causal questions blocked.

This is the motivating case study. AMY1/AMY2 dosage is **copy number**
(structural variation), not a SNP genotype, so all quantities below are integer
copy numbers. Every number here is reproduced by `analysis/amy.py`; commands
are given inline. Labels: `[observation]` = read directly from the supplied
tables; `[inference]` = my reading of them; `[speculation]` = explicitly not
established.

## 1. Tables, columns, individuals, copy-number distribution

`[observation]` `data/raw/amy/copy_number_table_for_supplement.tsv` is one row
per individual (5,667 data rows). Observed columns:

```
type  source  sample  p1  p2  AMY1  AMY2A  AMY2B
```

- `type` in {Modern, Ancient, Archaic, GTEx, Trio}
- `source` = publication / cohort provenance (e.g. `1KG`, `HGDP`,
  `Allentoft_Sikora_etal2024`, `Marchi_etal2022`, `Margaryan_etal2020`)
- `sample` = individual identifier
- `p1`, `p2` = population / region labels (**observed** metadata strings, not
  model output)
- `AMY1`, `AMY2A`, `AMY2B` = integer copy numbers from the Sudmant-lab
  graph-genotyping pipeline. These are **observed calls**, not inferred by this
  module.

`[observation]` Cohort sizes and AMY1 distribution (from `python3 analysis/amy.py --describe`):

| type | n | AMY1 mean | AMY1 range |
|---|---|---|---|
| Modern | 3,594 | 7.35 | 2–20 |
| Ancient | 533 | 6.89 | 2–17 |
| Archaic | 4 | 2.00 | 2 (all) |
| GTEx | 838 | 6.19 | 2–18 |
| Trio | 698 | 7.24 | 2–20 |

`[observation]` N with ancient provenance = **537** (533 `Ancient` + 4
`Archaic`). Ancient sources: Margaryan_etal2020 (333), Allentoft_Sikora_etal2024
(151), Allentoft_etal2015 (35), Marchi_etal2022 (14). The four Archaic
individuals are Denisova + 3 Neanderthals (Altai, Chagyrskaya, Vindija), all
AMY1=2. `[inference]` The 533 modern-era "ancient" samples plus the archaic
outgroup are the only temporally ordered AMY observations available.

`[observation]` `AMY2A` ∈ {0..6} (mode 2), `AMY2B` ∈ {2..7} (mode 2). AMY1 is
the variable gene; AMY2 copy number is near-constant.

`[observation]` Dates/coordinates are **not** in the copy-number table. They
come from `ancient_metadata_for_supplement.tsv`
(`sample, p2_l0, p2_l1, site, country, region, age, latitude, longitude,
data_source`; 533 ancient rows, `age` in years BP, range 250–11,627.5 BP,
median 1,075 BP). Structural haplotype labels come from
`haplotype_structure_for_supplement.tsv`
(`type, source, sample, p1, p2, best_matching_score, H1, H2`; 4,580 rows, but
only **288** ancient rows — a coverage gap vs. the 537 ancient copy-number rows).
Both are read directly out of `amylase_diversity_project-v1.1.zip` by the loader.

## 2. Linkability of ancient AMY individuals to Akbari / AADR

Command (reproduces every number below):

```
python3 analysis/amy.py --overlap
```

`[observation]` Results:

- `n_ancient_amy`: **537**
- `n_aadr_matched`: **522** — ancient AMY samples matched to an AADR v66.p1
  individual (by genetic id, individual id, or skeletal code)
- `n_aadr_matched_with_date_coord`: **522** — all AADR-matched ancients carry a
  date mean (BP) and latitude/longitude. `[inference]` So **522/537 = 97.2%**
  of ancient AMY individuals are dateable and geo-locatable via AADR.
- `n_akbari_ena_runs`: 10,016 sample identities parsed from
  `ENA_PRJEB106907_runs.tsv` (`submitted_ftp` bam names)
- `n_akbari_linked_to_aadr`: 384 — Akbari ENA identities that resolve to an AADR
  v66.p1 row
- `n_amy_ancient_in_akbari`: **0**

`[observation]` The **overlap for the trajectory test is 0**: none of the 537
ancient AMY individuals can be linked to an Akbari individual from the tables we
hold. `[inference]` Two independent naming/coverage gaps cause this:

1. Ancient AMY samples use Allentoft_Sikora2024 `NEO*` ids and older
   `Margaryan`/`Marchi` ids; Akbari identities are `I*` ids (Reich-lab style).
2. Akbari's set is dominated by 9,568 newly-reported individuals; only 384 of
   its 10,016 identities resolve to AADR v66.p1 at all, and those 384 are
   disjoint from the AMY ancients.

`[inference]` Consequence: a trajectory-label-based AMY test (using Akbari
trajectory labels) **cannot be run from the supplied data**. A date-only AMY
trend over the 522 AADR-linked ancients *is* possible and is supported by
`amy1_trend`, but it is a single pooled time series, not per-trajectory.

## 3. The five project questions

Answers use the trajectory table contract from `analysis/trajectory_infer.py`
(`individual, date, trajectory_id`). **That file does not exist yet**, so all
replication questions are blocked: `analysis/amy.py` refuses to invent labels
and exits non-zero (verified).

1. **Does AMY copy number rise along multiple independent trajectories?**
   *TBD — blocked on `analysis/trajectory_infer.py`* (no trajectory labels
   exist). Descriptive prerequisite satisfied: 522/537 ancient AMY individuals
   have dates and coordinates, so a per-trajectory trend is computable *once
   labels exist*.
2. **Are rises confined to replacements?** *TBD — blocked on trajectory labels*
   (needs `trajectory_id` plus an ancestry/replacement covariate, neither
   present).
3. **Do several structural haplotypes rise?** *TBD — blocked on trajectory
   labels*. Descriptive prerequisite: `[observation]` 288 ancient individuals
   have structural haplotype labels (`H1`/`H2`), spanning **10 distinct
   structural haplotypes** (dominant `H3^r` 215, `H5` 115, `H1^a` 96, `H7` 56 …).
   Whether several of these rise over time cannot be answered without
   trajectory/time-resolved labels; a haplotype-frequency-vs-date trend is
   computable for the 288 labeled ancients but is not a replication test.
4. **Do independent mutational backgrounds converge?** *TBD — blocked on
   both trajectory labels and sequence-level mutational-background data*
   (the archives here carry copy number and structural haplotype, not the
   flanking SNP haplotypes needed to define independent backgrounds).
5. **Does the evidence survive excluding close trajectories?** *TBD — blocked
   on trajectory labels and on a relatedness/trajectory-close pair table*
   (no such table is supplied).

## Reproduction

```
python3 analysis/amy.py --describe   # N individuals + copy-number distributions
python3 analysis/amy.py --overlap    # all counts in section 2 (overlap = 0)
python3 analysis/amy.py --trend TRAJ.tsv   # AMY1 slope per trajectory_id
python3 analysis/amy.py --figures    # fig/amy1_copy_number_distribution.png, fig/amy_overlap.png
```

`--trend` exits 2 with a loud message if the trajectory table is missing.

## Independent anchored-deme CN audit — appended 2026-10-08

[observation] New analysis/directional_consistency.py uses this module's CN loader and AADR biological identity linkage rather than requiring ENA/Akbari overlap. **522** ancient rows have actual numeric AADR date/latitude/longitude; **424** survive genotype/identity rebuilding criteria, **87** enter **13** qualifying long AADR paths. Only **one** has >=3 CN time nodes with >=5 carriers/node across >=2,000 years. The prespecified between-successive-demes feasibility check finds **zero** >=3-deme nonoverlapping successions with >=5 carriers/deme in one geographic cell.

[inference] The earlier zero Akbari-ENA overlap does NOT block AADR trajectory labeling; the older overlap helper also used truthiness of wrong metadata columns, which the new analysis validates against columns 10/17/18. Insufficient CN longitudinal replication, not Akbari naming, is the measured limitation now.

[observation] Held-out SNP control pass fraction **13/208=6.25%** (95% interval **3.37–10.45%**) passes the preliminary process gate. **CN null remains uncalibrated**; no identifiable drift-plus-jump CN ordering probability is reported. See docs/directional-consistency-results.md for all CN node means/slopes and limitations. [inference] These data cannot prove or disprove the accepted parallel-AMY argument; SNP Ne is not a CN mutation/haplotype model.

---

# CORRECTED AMY ANALYSIS (appended 2026-10-08)

**This supersedes the "TBD - blocked on trajectory_infer.py" conclusions above.** The
error was ours: we analysed AMY from `copy_number_table_for_supplement.tsv` plus our own
AADR-mapped ancient subset, and treated the paper's ancestry labels as mere annotation.
The paper's own supplement already contains the structure-aware ancient analysis, and we
never opened the released analysis archive (`amylase_diversity_project-v1.1.zip`) that
was sitting in `data/raw/amy/`.

## The paper's own data answers the question directly

`graph_genotyping/assess_graph_genotypes/supplementary_tables/ancient_metadata_for_supplement.tsv`
has **533 ancient individuals** with `p2` ancestry components, dates, and coordinates.
Joining to the CN table (533/533 join):

| ancestry (p2) | n | mean AMY1 | mean age (BP) |
|---|---:|---:|---:|
| Iron age to early modern | 333 | 7.25 | 934 |
| Bronze age | 46 | 6.63 | 3,893 |
| Steppe pastoralists | 16 | **8.25** | 4,068 |
| Neolithic farmers | 34 | **7.15** | 5,121 |
| EHG | 49 | 5.39 | 6,861 |
| WHG | 32 | 5.75 | 7,694 |
| Early farmers | 22 | **5.45** | 8,414 |
| CHG | 1 | 8.00 | 9,724 |

**[observation] The direction is exactly as the original argument claimed**: the older
components (WHG 5.75, EHG 5.39, Early farmers 5.45 - all around 5.4-5.8) are lower than
the later ones (Neolithic farmers 7.15, Steppe pastoralists 8.25, Bronze age 6.63, Iron
age+ 7.25). This is computable straight from the released supplement.

Global replication is also present in the same release
(`hdgp_tgp_AMY_popgen/CN_metadata_filtered_withSubsistence.txt`, 3,395 individuals):
AMY1 by genetic region is highest in **OCE 8.72, EAS 8.08, AMR 7.58** and lowest in
**EUR 6.89, AFR 6.71**; and within the Americas, **PEL (Andean) 9.52 and Pima 10.29**
are the highest populations in the dataset versus **Surui (Amazonian) 5.25**.

## But the naive "independent innovation" claim is NOT supported by the haplotype labels

`haplotype_structure_for_supplement.tsv` gives per-individual structural haplotype labels,
and CN is strongly haplotype-dependent:

| haplotype | carriers | mean AMY1 |
|---|---:|---:|
| H9 | 104 | 13.95 |
| H7 | 529 | 11.16 |
| H5 | 2,414 | 9.01 |
| H4A2B2 | 278 | 8.27 |
| H3^r | 4,127 | 6.49 |
| H2A2B2 | 141 | 5.43 |
| H1^a | 949 | 4.30 |

**[observation] All 10 ancient haplotype labels also occur in moderns** (only `H3A3B3` is
modern-only). So high copy number tracks *specific shared structural haplotypes*, and
their frequencies differ sharply between populations - H5 is 0.44-0.50 in CHS/CHB/JPT but
only 0.16-0.19 in TSI/IBS/YRI. Haplotype-frequency differences are therefore a plausible
**carriage** route for the elevated Andean/East-Asian copy number, which is precisely the
confound the Akbari-lab objection raises.

Ancient composition moves in the same direction as CN: WHG is dominated by low-CN
`H2A0`/`H1^a`, while later groups carry more `H5`/`H7`/`H4A2B2`.

## Honest status

- **Corrected:** "AMY is not testable" was wrong. The paper's own supplement contains the
  structure-aware ancient analysis, it is on disk, and it shows the increase.
- **Not resolved by labels alone:** whether the high-CN haplotypes in the Americas are
  *independently derived* or *carried* from a shared ancestral pool. Shared labels are
  consistent with carriage; the paper's "recurrent evolution" claim must therefore rest on
  the finer structural/sequence-level analyses in the archive (bundle/structure trees,
  `haplotype_all_structures*`, PGGB graphs), which we have NOT yet examined.
- **The decisive test remains donor-aware**: given a steppe/EHG donor whose haplotype
  composition already has higher mean CN, how much of the later West Eurasian increase is
  predicted by admixture proportion alone versus requiring within-component change? That is
  answerable with `p2` proportions plus the haplotype-specific CN table, and it is the test
  the argument actually needs.

## CN-as-SNP counterfactual — appended 2026-10-09 (new analysis, prior conclusions retained)

[observation] `docs/cn-as-snp-prereg.md` preceded fitting; `docs/cn-as-snp-results.md` contains the complete same-cohort controls, five perfect-tag haplotype dosages, continuous AMY1, GRM-on/off contrasts and 199-attempt label-permutation nulls. Unlike the earlier ENA-only linkage result, the released imputed bundle directly matches **225/288** labelled ancients by unique biological-ID prefix. This is a different data source/denominator, not a rewrite of the historical overlap claim.

[observation] The supervisor-approved ordering gate passed (GRM-on |Z| LCT=5.918, SLC45A2=9.297, released lead AMY SNP=0.676). H1^a |Z|=4.116 on / 4.246 off and H5=2.706 both; no GRM-collapse occurs with the 1,133-marker genome-wide LOCO construction. Neither clears the directly transported 3.16 X-equivalent cutoff. The uncertain same-cohort LCT-extrapolated |X| values are 11.493 and 7.556, respectively; these are not reproduced Akbari statistics. All five ancestry-conditioned haplotype P values exceed 0.05, while unrestricted P is 0.005 for H1^a and 0.010 for H5. H9's primary GRM-off fit failed, and a separately labelled longer-iteration same-model diagnostic is reported.

[inference] Perfect tags strengthen the trend relative to the released lead SNP, but the original-scan counterfactual remains unresolved: transported calibration is not faithful, ancestry-associated change remains compatible, late ancients dominate and the sparse GRM may not represent the original relatedness correction. This does not establish over-correction, independent repeated selection or a real-scan copy-number call.

## NEW LOCAL-ANCESTRY STAGE 1 — 2026-10-09T17:03Z

[observation] See `docs/local-ancestry-prereg.md` and `docs/local-ancestry-results.md`. Implemented an unphased 15-state ancestry-pair HMM with 22,986 AIMs and nine sensitivity settings. Painted 278 structural-haplotype targets (239 with supported AMY vicinity), but validation FAILED: neither IRAN holdout is adequately covered; EHG median correct-source fraction=.514<.60; window/penalty proportions are unstable. Covered-reference dominant classification19/19 and negative-control WHG=.958 do not waive these failures. AMY local-conditioned regression and outcome frequencies were STOPPED, not reported.

[inference] No validated answer on fixed-ancestry AMY persistence is available; prior carriage-dominated conclusion is neither changed nor corroborated. Raw CN means and released AMY2204/max|X|3.794956889 were rechecked; previous225-person global-label slope=.048575 versus floor=.145764 and carriage94.61% reproduce. New strict unique-biological AADR CN mapping gives519 dated-coordinate matches versus supplied522, explicitly reported as a three-person disagreement. A future revision needs separately preregistered coverage/segment validation and phase-linked structural alleles, not retroactive tuning.

## APPENDED full-imputed-panel haplotype backgrounds — 2026-10-09

[observation] New preregistered full-panel lane uses ALL16,339 imputed individuals (32,678 unphased chromosome slots), all eligible window SNPs, and all9160 supplied structural H1/H2 slots. Imputed usable coverage is100% in both windows at all five loci; no curated ancestry references or missing-source holdouts. Annotation auditing found the initial AMY exclusion missed AMY1 copies through104301314; supervisor-approved post-results correction generalises transcript-complete exclusions to every locus, with all initial outputs retained as INVALID-exclusion diagnostics.

[inference] Corrected backgrounds still fail known-structure and stability gates; AMY and canonical conditioned interpretations STOP. Full-panel coverage fixes the old coverage bottleneck but is not evidence of validated backgrounds, phase or power. Input structural identity inventory and outcome-blind design floors are reported without claiming molecular recurrence or a validated within-background AMY effect. See `docs/haplotype-background-prereg.md`, `docs/haplotype-background-results.md` and `data/derived/haplotype_background/`; prior authors' conclusions are unchanged.
