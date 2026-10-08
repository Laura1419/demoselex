# AADR v66.p1 QC report

Scope: quality control of the AADR v66.p1 1240K (`v66.p1_1240K.aadr.patch.PUB`)
and Human Origins (`v66.p1_HO.aadr.patch.PUB`) bundles, as input to a downstream
trajectory analysis. All numbers below are reproducible from the commands and
seed recorded here and are mirrored in `data/derived/qc_manifest.json`.

## Provenance / reproducibility

- Commands: `python3 -m analysis.load_aadr` (self-test) and
  `python3 analysis/qc.py` (writes `data/derived/qc_manifest.json`, ~1m32s).
- Seed (site and individual sampling, relatedness): **20261008**.
- Versions: Python 3.12.3, NumPy 1.26.4 (also stored in the manifest).
- Input files are read-only and unmodified.

## File format (verified)

`[observation]` The `.geno` files are not bare EIGENSTRAT. They begin with a
48-byte `TGENO` header:

```
TGENO   23089 1233013 5b63c3d2 5ccb0c38   # main  (NUL-padded to 48 bytes)
TGENO   27594  584131 8c17d6d1 80974215   # HO
```

The payload is **individual-major**, 2 bits/site, 4 sites per byte, one
contiguous row of `ceil(n_snp/4)` bytes per `.ind` line. This is confirmed by the
exact byte-size identity (`[observation]`):

| bundle | n_ind | n_snp | 48 + n_ind*ceil(n_snp/4) | on-disk size | match |
|---|---|---|---|---|---|
| 1240K patch PUB | 23,089 | 1,233,013 | 7,117,276,654 | 7,117,276,654 | yes |
| HO patch PUB | 27,594 | 584,131 | 4,029,634,650 | 4,029,634,650 | yes |

Genotype codes are 0/1/2 for the diploid state in `.snp` ref/alt orientation and
3 = missing. `[inference]` The ref/alt orientation is taken from `.snp` columns
5/6 as documented for EIGENSTRAT; we did not independently verify that AADR's
A1/A2 allele ordering matches, and MAF below is therefore defined on the
`snp_alt` orientation value.

`.ind` = 3 columns `(id, sex, group)`; `.snp` = 6 columns
`(rsid, chrom, genetic_pos, phys_pos, ref, alt)`; header counts equal the file
line counts (`[observation]`, self-test asserts this).

## Counts

| | value |
|---|---|
| Individuals (1240K patch, `.ind` rows) | 23,089 |
| Sites (1240K patch, `.snp` rows) | 1,233,013 |
| Individuals (HO subset) | 27,594 |
| Sites (HO subset) | 584,131 |
| `.anno` rows (main) | 23,089 |

## Date semantics and coverage (`v66.p1_1240K.aadr.PUB.anno`)

`[observation]` The AADR ancestor file documents the primary time axis as
**"Date mean in BP in years before 1950 CE"** (OxCal mu for a direct
radiocarbon date; average of range for a contextual date), with a companion
standard deviation in BP and a free-text "Full Date" column that mixes formats
(`calBCE`, `BP`, `BCE`). **Convention adopted:** use the numeric
`Date mean in BP` column directly as "years before 1950 CE"; do not re-calibrate
or convert the `Full Date` strings, and **flag** rather than guess every row
whose method is `Modern`, whose method text carries `WARNING`/`CAUTION`, or whose
method field is anomalous.

| quantity | value |
|---|---|
| rows with numeric `Date mean BP` | 23,089 (all rows) |
| rows with numeric latitude/longitude | 22,305 |
| rows missing coordinates | 784 |
| rows with `Modern` date method | 3,970 |
| rows with usable date **and** coordinates **and** not `Modern` | **19,029** |
| rows whose date method carries WARNING/CAUTION | 1,647 |
| date BP min / median / max (usable) | 60 / 2,554 / 185,000 |
| date BP quantiles p5/p25/p75/p95 (usable) | 593 / 1,225 / 4,767 / 7,725 |

Method categories `[observation]`: Context 9,482; Direct 9,175; Modern 3,970;
Other/Unlabeled 462 (a handful of Jinlianshan rows carry a date *range* string in
the method column and are treated as unlabeled).

`[observation]` Outliers that must NOT be treated as modern-human samples:
11 rows have `date_bp > 50,000`, of which the extremes are Neanderthal/Denisovan
references, e.g. `Denisova25.SG` (185,000 BP, "Context: molecular dating"),
`AltaiNeanderthal.DG` (110,450), `Denisova11.SG` (98,700). At the other end,
`Khwit.SG` has `date_bp = -4` (a 20th-century known individual). These are
flagged by the policy below.

`[inference]` Because "Context" dates are ranges averaged to a mean and
"Direct" dates are calibrated densities, the `date_sd_bp` column is not
homogeneous across rows; a downstream model should treat date uncertainty as
method-dependent, and rows with WARNING/CAUTION or unlabeled methods as
lower-confidence.

## Missingness and MAF (main 1240K set)

`[observation]` Computed on a fixed random subset of **200,000 sites** (seed
20261008) in one streaming pass over the memory-mapped matrix. Individual
missingness = fraction of sampled sites coded 3.

| quantity | p0 | p25 | p50 | p75 | p95 | p100 |
|---|---|---|---|---|---|---|
| per-individual missingness | 0.0004 | 0.163 | 0.422 | 0.810 | 0.983 | 1.000 |
| per-site missingness | 0.000 | 0.382 | 0.452 | 0.556 | 0.693 | 1.000 |
| per-site MAF | 0.000 | 0.048 | 0.176 | 0.329 | 0.465 | 0.500 |

- mean per-individual missingness = 0.478.
- individuals with missingness > 0.90: **3,997**.
- among the 19,029 date+coord-usable individuals, missingness median = 0.535,
  p95 = 0.986 `[observation]` — i.e. many date-usable samples are low-coverage.
- monomorphic sites in the subset: 281.
- sites with MAF >= 0.05: **149,123 / 200,000 (74.6%)**.

`[inference]` The high median site missingness (~45%) is expected: the 1240K
"patch" matrix legitimately contains `3` (missing) at sites each sample does not
cover, and AADR v66 mixes high-coverage genomes with many low-coverage
captures.

## Relatedness screen (HO subset)

`[observation]` Crude IBS screen on the HO bundle: sampled **1,200 individuals**
and **15,000 sites** (seed 20261008) and computed, per pair, the proportion of
jointly non-missing sites with the same genotype ("IBS concordance") and the
proportion of opposite-homozygote sites ("IBS0"). Pairs sharing fewer than
**1,000** jointly-called sites were dropped (186,707 of 719,400 pairs).

| quantity | value |
|---|---|
| pairs tested (>=1,000 shared sites) | 532,693 |
| IBS concordance p25 / p50 / p75 / p95 | 0.623 / 0.647 / 0.756 / 0.772 |
| IBS0 p25 / p50 / p75 / p95 | 0.126 / 0.153 / 0.237 / 0.255 |
| pairs with IBS >= 0.90 (related candidate) | **1** |
| pairs with IBS >= 0.97 (duplicate candidate) | **1** |
| pairs with IBS0 < 0.005 (near-identical) | **1** |

`[observation]` Exactly one pair is simultaneously the IBS>=0.97 and the
IBS0<0.005 pair — a clear duplicate / same-individual pair (or the same genome
entered under two `.snp`-suffix IDs). IBS0 ~ 0 while IBS ~ 1 is the signature of
a duplicate, not of first-degree relatives.

`[speculation]` The screen is deliberately crude and subsampled; it is powered
to catch duplicates and close pairs, not to estimate a kinship matrix. Regions
of isolation / strong structure can inflate IBS between unrelated individuals,
so the 0.90 threshold is conservative and the absence of pairs just below it
should not be read as "all unrelated".

## Exclusion policy for a downstream trajectory analysis

`[observation]`-driven policy; downstream code should drop, and log why:

**Individuals — drop if any of:**
1. No numeric `Date mean BP`, no latitude, or no longitude (784 lack coords).
2. `date_method` starts with `Modern` (3,970) — not part of an ancient
   time-series.
3. `date_bp` outside a modern-human window; concretely drop the 11 archaic /
   reference samples with `date_bp > 50,000` and the one negative
   (`Khwit.SG`, -4 BP).
4. Per-individual missingness > 0.90 (3,997 individuals) — genotype calls too
   sparse to place reliably.
5. One member of every duplicate/related pair above the IBS threshold (>= 1
   pair at >= 0.97 on the HO screen). Keep the higher-coverage / better-dated
   member.
6. `[inference]` Rows with `ASSESSMENT` not `Pass`, or method WARNING/CAUTION
   (1,647), should be flagged as low-confidence; drop them unless the analysis
   explicitly models date uncertainty.

Only 19,029 individuals pass filters 1–2 up front; the intersection with the
missingness filter (3) and the relatedness filter (5) is not persisted in this
pass and should be computed once the per-individual missingness vector is
written out (see residual gap below).

**Sites — drop if any of:**
1. Monomorphic in the analysis cohort (subset estimate: 281/200,000 ~ 0.14%,
   but recompute on the retained cohort, not the full 23,089).
2. MAF < 0.05 (subset estimate: 50,877/200,000 = 25.4% drop) — keep the 74.6%
   with MAF >= 0.05 for a common-variant trajectory.
3. `[inference]` Sites with unusually high missingness (> cohort median +
   margin, or > 0.7) to reduce imputation burden; the subset p95 is 0.693.

Use the HO bundle only when a variant-overlap screen is needed; it is a
different (larger, lower-coverage) individual set and must be merged on `id`,
not assumed row-aligned with the 1240K set.

## Residual gaps / risks

- `[observation]` The single duplicate/related pair's identity is not persisted
  in the manifest (counts only); rerun with ID logging before applying filter 5.
- The date `date_bp` field is taken at face value; WARNING rows (1,647) and the
  7 Jinlianshan range-in-method rows were flagged but not re-parsed.
- Relatedness is a subsample screen (1,200 x 15,000), not a full kinship
  estimate; it will miss weaker relatives.
- Per-individual missingness vector is not written to disk, so the
  `date-usable AND coverage-usable` intersection remains to be computed.
