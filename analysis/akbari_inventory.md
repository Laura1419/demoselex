# Akbari et al. 2026 selection statistics — data inventory

> Lane: data loading + inventory only. GLMM reimplementation lives in `sim/akbari_glmm.py`.
> Loader: `analysis/load_akbari.py`. Derived tables: `data/derived/akbari_stats.json`,
> `data/derived/akbari_signal_absXge4.tsv`, `data/derived/akbari_background200k.tsv`,
> `data/derived/akbari_top.tsv`.
>
> Source file: `data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz`
> (674,584,607 B, sha256 prefix `b888d5a4`, CC0, Akbari et al. 2026 *Nature* 654:419-428,
> doi:10.1038/s41586-026-10358-1). Environment: python 3.12.3, pandas 2.1.4, numpy 1.26.4.
> All numbers below are reproducible from the loader + the exact commands in `## Reproduce`.

## Row count and chromosome coverage — [observation]

- **9,739,624 data rows** (streamed count). Matches the paper's stated final variant set
  (9,739,624 variants) exactly; the loader asserts this in `--inspect`.
- Autosomes **1–22 only**, all present. Per-chromosome counts range from **124,043** (chr21)
  to **818,116** (chr2). No chrX/Y/MT rows.
- Build: **GRCh37 / hg19** (stated in the in-file header).
- `CHROM` and `POS` are the only positional keys; `ID = CHROM_POS_REF_ALT`.

## Columns: verified meaning and units — [observation]

Header lines are shipped *inside* the gzip as `##COLUMN=<Description=...>`; the descriptions
below are quoted verbatim from them (authoritative source of units). Types are the dtypes
pandas infers from a real chunked read.

| # | Column | dtype | Meaning (verbatim from file header) | Units |
|---|--------|-------|--------------------------------------|-------|
| 1 | `CHROM` | Int64 | Chromosome (GRCh37/hg19) | 1..22 |
| 2 | `POS` | int64 | Position on the chromosome (GRCh37/hg19) | bp, 1-based |
| 3 | `REF` | str | Reference allele (GRCh37/hg19) | nt string |
| 4 | `ALT` | str | Alternative allele | nt string |
| 5 | `ANC` | str | Ancestral allele | nt string; `NA` when unknown |
| 6 | `ID` | str | Unique identifier = `CHROM_POS_REF_ALT` | string |
| 7 | `RSID` | str | Reference SNP ID | string; `NA` when absent |
| 8 | `AF` | float64 | Allele frequency of the alternative allele (ALT) | fraction in [0,1] |
| 9 | `S` | float64 | Selection coefficient of the alternative allele (ALT) | see note |
| 10 | `SE` | float64 | Standard error for the selection coefficient | same as `S` |
| 11 | `X` | float64 | Selection X statistic (Z-score adjusted for inflation) | dimensionless Z |
| 12 | `P_X` | float64 | P-value of X, assuming normal | probability |
| 13 | `POSTERIOR` | float64 | Posterior probability of a true selection signal given X | probability |
| 14 | `FDR` | float64 | Expected proportion of false positives among detected signals given X | probability |
| 15 | `CHI2_BE` | float64 | Chi-square for batch effect (higher = more suspect) | chi-square |
| 16 | `FILTER` | str | Status / warning / error messages | see below |

**`FILTER` values actually present** [observation]:
`PASS` = 8,804,884 rows; `WARNING:MAF<0.01` = 934,740 rows. No error-level flags appear.

### Which column is the time-slope / reported Z statistic? — [observation / inference]

- **`X` is the statistic the paper reports as genome-wide significant.** The header defines it
  as *"Selection X statistic (Z-score adjusted for inflation)"*; it is the inflation-corrected
  Z of the per-variant time slope. The paper's genome-wide threshold is stated as
  **|Z| > 5.45** (see `docs/C-akbari-reconstruction.md`, from the GWAS-enrichment calibration),
  and that threshold applies to `X`. `P_X` is the two-sided normal p-value of `X`
  (verified: for `X = 1.6698`, `P_X = 0.0950 = 2*(1-Φ(1.6698))`; for `X = -16.22`,
  `P_X = 3.4e-59`).
- **`S` is the fitted time slope itself** (the "selection coefficient of the ALT allele"),
  and `SE` is its standard error. `X` is the corresponding inflation-adjusted Z.

### Ambiguity that remains — [observation]

- The file does **not** state the units of `S`/`SE`. In the paper the fixed effect is a
  **time slope** (`β_time`) that "captures correlation between allele frequency and time after
  correcting for population structure" and "differs from classical *s* in discrete-generation
  population-genetics models" (`docs/C-akbari-reconstruction.md`). The observed magnitudes
  (|S| ~ 1e-3, e.g. 0.0019, 0.0179) are consistent with a per-generation frequency slope, but
  the exact time unit (per generation vs per year vs per kyr) is **not verified from the file**.
  Do not treat `S` as a classical selection coefficient without checking Supplementary
  Table / Methods. `X`'s units (dimensionless Z) are unambiguous.
- `X` is **not mean-zero**: mean = **+0.412**, std = 1.491, **62.3% of variants have X > 0**,
  37.7% have X < 0. [observation] This is the genome-wide distribution of the *signed* statistic
  for the ALT allele. The positive shift is large relative to a null N(0,1) and is **not**
  explained by the "inflation adjustment" wording alone. [speculation] It may reflect (a)
  genuine genome-wide directional allele-frequency change over the sampled time span, (b) a
  sign-convention effect (ALT = non-reference in hg19, often the derived allele, plus
  ascertainment to polymorphic sites), or (c) an uncentred statistic definition. This needs a
  look at the Methods before `X` is used as a two-sided Z in simulations — flagging for the
  method lane rather than resolving here.

## Statistic distribution (`X`), n = 9,739,624 — [observation]

Exact (full-column) values; seedless, deterministic.

| stat | value |
|------|-------|
| min / max | -16.224 / +16.643 |
| mean / std | +0.4120 / 1.4905 |
| frac > 0 / < 0 | 0.6232 / 0.3768 |
| `X` quantiles 0.1/1/5/25/50/75/95/99/99.9/99.99 % | -6.070, -3.321, -1.985, -0.504, **+0.429**, +1.358, +2.742, +3.878, +6.258, +11.100 |
| `|X|` quantiles 50/90/99/99.9/99.99/99.999/99.9999 % | 0.981, 2.470, 4.267, 7.243, 11.985, 15.577, 16.111 |
| # variants `|X| >= 3` | 481,989 |
| # variants `|X| >= 4` | **131,103** |
| # variants `|X| >= 5` | 46,381 |
| # variants `|X| >= 5.45` (paper threshold) | **31,939** |
| # variants `|X| >= 6` | 21,988 |
| # variants `|X| >= 7` | 11,590 |
| # variants `|X| >= 8` | 6,180 |
| # variants `|X| >= 8` … max | 6,180 (max +16.64) |

### Genome-wide-significance threshold actually present in the data — [observation]

- The `FDR` and `POSTERIOR` columns are **floored**: `FDR` never drops below **0.0099** and
  `POSTERIOR` never exceeds **0.9901** for any row. So no row carries an FDR of exactly 0.
- Using the data's own FDR column, **`FDR <= 0.05` first appears at `|X| ~ 5.02`** (min `|X|`
  among sampled rows with `FDR <= 0.05` was 5.0228); by `|X| >= 5` the **median FDR is already
  0.0099 (the floor) and median POSTERIOR is 0.9901**. At `|X| >= 4` median FDR is 0.139 and
  median POSTERIOR 0.664.
- The paper's stated threshold is **|Z| > 5.45**, which here selects **31,939** variants.
  The data-derived FDR<=0.05 crossing (~5.0) and the paper's 5.45 bracket the same region;
  they are not identical, so record which one a downstream analysis intends to use.

## Top 25 loci by `|X|` — [observation]

Source: `data/derived/akbari_top.tsv` (written by `load_akbari.py --top 25`; the file and this
table are generated from the same tool).
"Nearest gene" is **not trivially available**: no gene/GTF/GFF annotation exists in the repo,
so no gene is asserted here. The clusters are given as `chrom:pos` bands. [speculation] The
chr2 band 136.60–136.71 Mb is in the well-known LCT/MCM6 region, and the chr5 (~131.6–131.8 Mb)
and chr11 (~61.5–61.7 Mb) loci recur throughout the very top of the ranking — but these are
geographic inferences from coordinates, not annotation lookups, and should be verified before use.

| Rank | variant id | chrom:pos | ALT freq | `X` | POSTERIOR | FDR |
|------|-----------|-----------|----------|-----|-----------|-----|
| 1 | 5_131664516_TTATATA_T | 5:131664516 | 0.3430 | 16.6433 | 0.990 | 0.0099 |
| 2 | 5_131640536_A_G | 5:131640536 | 0.4293 | 16.5297 | 0.990 | 0.0099 |
| 3 | 2_136608646_G_A | 2:136608646 | 0.2714 | 16.5232 | 0.990 | 0.0099 |
| 4 | 5_131797578_T_C | 5:131797578 | 0.5746 | -16.2242 | 0.990 | 0.0099 |
| 5 | 11_61603510_C_A | 11:61603510 | 0.4275 | -16.1534 | 0.990 | 0.0099 |
| 6 | 5_131710399_A_C | 5:131710399 | 0.6346 | -16.1431 | 0.990 | 0.0099 |
| 7 | 11_61588305_A_G | 11:61588305 | 0.4239 | -16.1411 | 0.990 | 0.0099 |
| 8 | 5_131707429_G_A | 5:131707429 | 0.6345 | -16.1237 | 0.990 | 0.0099 |
| 9 | 5_131709032_T_TAATA | 5:131709032 | 0.6345 | -16.1237 | 0.990 | 0.0099 |
| 10 | 11_61569354_A_AT | 11:61569354 | 0.4097 | -16.1144 | 0.990 | 0.0099 |
| 11 | 11_61593005_A_G | 11:61593005 | 0.4257 | -16.1093 | 0.990 | 0.0099 |
| 12 | 5_131705458_C_G | 5:131705458 | 0.6351 | -16.0994 | 0.990 | 0.0099 |
| 13 | 11_61596633_C_T | 11:61596633 | 0.4258 | -16.0988 | 0.990 | 0.0099 |
| 14 | 11_61600342_A_C | 11:61600342 | 0.5711 | 16.0740 | 0.990 | 0.0099 |
| 15 | 5_131723288_T_C | 5:131723288 | 0.4334 | 16.0654 | 0.990 | 0.0099 |
| 16 | 5_131723289_T_A | 5:131723289 | 0.4333 | 16.0564 | 0.990 | 0.0099 |
| 17 | 11_61592362_A_G | 11:61592362 | 0.4246 | -16.0541 | 0.990 | 0.0099 |
| 18 | 5_131677047_G_T | 5:131677047 | 0.3645 | 16.0446 | 0.990 | 0.0099 |
| 19 | 11_61604814_C_A | 11:61604814 | 0.4288 | -16.0306 | 0.990 | 0.0099 |
| 20 | 11_61609750_C_T | 11:61609750 | 0.4320 | -16.0190 | 0.990 | 0.0099 |
| 21 | 2_136616754_C_T | 2:136616754 | 0.2771 | 16.0169 | 0.990 | 0.0099 |
| 22 | 5_131672657_T_C | 5:131672657 | 0.3106 | 16.0034 | 0.990 | 0.0099 |
| 23 | 11_61573684_T_C | 11:61573684 | 0.4157 | -16.0030 | 0.990 | 0.0099 |
| 24 | 5_131630852_A_C | 5:131630852 | 0.3104 | 15.9973 | 0.990 | 0.0099 |
| 25 | 5_131676320_C_T | 5:131676320 | 0.3106 | 15.9903 | 0.990 | 0.0099 |

Alternatively signed-`X` maxima (top by `+X`, e.g. variants #1–3,14,15) and minima (top by
`-X`, e.g. #4–13) are the two tails; `--top n` uses `|X|` by default, `top(by='x')` gives the
signed top.

## Derived tables (bounded) — [observation]

Exact filters and row counts:

| file | filter / construction | data rows |
|------|-----------------------|-----------|
| `akbari_signal_absXge4.tsv` | **all** variants with `abs(X) >= 4`; all 16 columns kept | **131,103** |
| `akbari_background200k.tsv` | systematic random sample of the full table: exactly K=200,000 rows at global row positions `floor(start + step*i)`, i=0..K-1, with `step = N/K = 48.69812`, `start = rng(seed=42).uniform(0, step) = 37.690205`; columns ID, CHROM, POS, REF, ALT, AF, X, S, SE, FILTER | **200,000** |
| `akbari_top.tsv` | top 25 rows by `abs(X)` (`load_akbari.py --top 25`); ties broken by file order | **25** |
| `akbari_stats.json` | computed summary (row counts, per-chrom counts, threshold counts, exact `X` quantiles, filter counts, seed) | — |

The 131,103 `|X|>=4` rows are comfortably bounded (≈1.3% of the table); a stricter
`|X| >= 5.45` cut would give 31,939 rows if a smaller signal set is wanted.

## Caveats — [observation unless noted]

1. `FDR` and `POSTERIOR` are **floored** at 0.0099 / 0.9901 — they cannot express stronger
   significance, so they are not a continuous ranking; use `X` for ranking.
2. `X` is **not centred at zero** (mean +0.41, 62% positive) — see the ambiguity note above.
3. **Tied `X` values** exist: e.g. `5_131707429_G_A` and `5_131709032_T_TAATA` both have
   `X = -16.12374955270621`. Ties (identical float, and identical `AF = 0.634485…`) recur in the
   top ranks, so "top-N by X" is order-sensitive at the tail. [inference] This likely reflects
   the statistic being computed on a shared/windowed or imputed block, but it is not verified.
4. `AF` is the **ALT-allele** frequency (header), which in hg19 is the non-reference allele,
   not necessarily the derived allele. `ANC` is often `NA`.
5. No gene annotation in-repo → nearest-gene column is intentionally omitted (not guessed).
6. GRCh37/hg19 coordinates; do not mix with GRCh38 resources.
7. `WARNING:MAF<0.01` flags 934,740 low-frequency variants that are still present in the table;
   decide per-analysis whether to keep them (they are not excluded here).

## Reproduce — [observation]

```
python3 analysis/load_akbari.py --inspect
python3 analysis/load_akbari.py --top 25 --out data/derived/akbari_top.tsv
# full-table stats / derived signal + background tables were produced with the
# streaming scripts recorded in docs/I-akbari-reanalysis.md (seed 42 for the sample).
```
