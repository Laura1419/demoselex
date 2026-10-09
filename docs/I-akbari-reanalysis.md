# I. Akbari reanalysis

> Owner task: locus-akbari. Status: TODO.

## Four-quadrant classification
| Quadrant | Akbari | Rep  | Interpretation |
| --- | --- | --- | --- |
| Q1 | high | high | strongest consensus |
| Q2 | high | low  | single broad/demographic event |
| Q3 | low  | high | **candidate over-correction** |
| Q4 | low  | low  | weak/null |

## Low-Akbari / high-replication loci in depth

## Akbari summary-statistics inventory (data lane)

> Added by the `load-akbari` data lane. Loading + inventory only; no GLMM reimplementation
> (`sim/akbari_glmm.py` owns that). Full details and the top-25 table live in
> `analysis/akbari_inventory.md`; loader is `analysis/load_akbari.py`.

**File.** `data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz`, CC0, Akbari et al.
2026 *Nature* 654:419-428. Streamed row count = **9,739,624**, matching the paper's final
variant set exactly. GRCh37/hg19, autosomes 1–22 only, per-chrom counts 124,043 (chr21) to
818,116 (chr2). Env: python 3.12.3 / pandas 2.1.4 / numpy 1.26.4.

**Reported statistic [observation].** The per-variant test statistic is `X` — header:
*"Selection X statistic (Z-score adjusted for inflation)"*. `S` is the fitted time slope
(selection coefficient of the ALT allele) and `SE` its standard error; `X` is that slope's
inflation-adjusted Z. `P_X` is the two-sided normal p-value of `X` (spot-checked against Φ).
The paper's genome-wide cut **|Z| > 5.45** applies to `X`. **Residual ambiguity:** the file does
not state units for `S`/`SE`; the paper calls `β_time` a time slope that "differs from classical
*s*". Do not treat `S` as a classical selection coefficient without the Methods / Supp. Table 4.
Only `X` is unambiguous as a dimensionless Z.

**Distribution [observation].** `X` is **not mean-zero**: mean +0.412, std 1.491, 62.3% of
variants positive, median +0.429, min/max -16.224/+16.643. Threshold counts: `|X|>=4` = 131,103;
`>=5` = 46,381; `>=5.45` = 31,939; `>=6` = 21,988. `FDR` is floored at 0.0099 and `POSTERIOR`
at 0.9901 (no row reaches FDR 0). Using the file's own FDR column, `FDR<=0.05` first appears at
`|X| ~ 5.02`, i.e. the data-derived cut and the paper's 5.45 bracket the same region but are not
identical.

**Derived tables [observation]:** `data/derived/akbari_signal_absXge4.tsv` (all 131,103 rows
with `|X|>=4`), `data/derived/akbari_background200k.tsv` (systematic random sample, K=200,000,
seed 42, step 48.69812), `data/derived/akbari_top.tsv` (top 25 by `|X|`),
`data/derived/akbari_stats.json` (summary).

**Caveats.** (1) FDR/POSTERIOR floored -> not a continuous ranking; rank by `X`.
(2) `X` centred at +0.41, not 0 — [speculation] genuine directional shift and/or ALT-ascertainment
sign convention; needs Methods before use as a two-sided Z. (3) Tied `X` floats at the tail
(e.g. two chr5 variants both `-16.12374955270621`) make "top-N" order-sensitive. (4) `AF` is the
ALT (hg19 non-reference) frequency; `ANC` often `NA`. (5) No in-repo gene annotation, so no
nearest-gene column is asserted. (6) 934,740 rows carry `WARNING:MAF<0.01` but are retained.

**Reproduce.**
```
python3 analysis/load_akbari.py --inspect
python3 analysis/load_akbari.py --top 25 --out data/derived/akbari_top.tsv
python3 -c "import analysis.load_akbari as L; print(len(L.load(chrom=22, min_abs_x=6)))"  # 177
```


## Third-pass reference-ranking comparison — 2026-10-09 (append only)

[observation] New bounded ranking uses uncentred Akbari |X| as a reference, never ground truth. Of108 scored variants with |X|>5.45,73=67.59% enter our top decile; the cache contains31,939 such variants, so this is NOT genome-wide recall. Whole-panel Spearman=-.1182. The panel is Akbari-ascertained (all |X|>=3.16 plus sampled lower-X cache) and percentiles refer only to its10,641 scored variants.

[observation] There are992 top-decile candidates NOT Akbari-significant; all identities/components are in data/derived/prioritization/top_decile_not_akbari_significant.tsv. Examples include rank20 10_113634916_A_AC (|X|1.636) and rank27 18_59876495_A_G (|X|.671). Detailed counts, shape caveats and prioritisation interpretation: docs/prioritization-results.md.

[inference] These are extra downstream phenotype-work candidates, NOT proven Akbari false negatives or evidence that structure correction suppressed true selection. The earlier hypothesis tests above remain intact and answer a different question.
