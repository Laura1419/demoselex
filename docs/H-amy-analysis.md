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
