# Results — between-history ordered rise (traj-trend)

> Owner lane: traj-trend. Producer: `analysis/traj_trend.py`.
> Pre-registration: `docs/traj-trend-prereg.md` (written before the run).
> Artifacts: `data/derived/traj_trend.tsv` (per-locus), `fig/traj_trend.png`.
> Labels: `[observation]` measured here; `[inference]` conclusion from it;
> `[speculation]` unverified.

## What was run (exactly)

`[observation]` `.venv/bin/python analysis/traj_trend.py --max-variants 8000`
completed in ~9 s and wrote `data/derived/traj_trend.tsv` incrementally (8 chunks
of 1000). Panel analysed: **8,000 variants = 2,000 candidates (`|X| >= 3.16`) +
6,000 background (`|X| < 3.16`)**, drawn from `data/derived/tier2_pool_cache.pkl`
(538,274 QC-confirmed Akbari variants; copied from `/tmp` so it survives cleanup).
Units: **m = 6 merged replication groups** (`merged_group` G0..G5; `N_eff =
m^2/sum(R) = 5.086`) with the matching 6x6 `traj_relatedness_v3.tsv` matrix; 5,737
of the trajectory individuals are present in the Akbari pgen bundle. The statistic
is the pre-registered relatedness-aware GLS slope of history mean allele frequency
on history placement `t = -mean(date_BP)` (larger `t` = more recent), per
`docs/traj-trend-prereg.md` §3.

Validation run `--max-variants 2000` also completes and appends rows (500
candidates + 1,500 background).

## Pre-registered decision rule (restated)

Primary empirical P = matched-control empirical P (background loci matched on
chromosome, Akbari AF, missingness, >= 1 Mb apart; prereg §4b), the calibrated
null; the fixed-scale permutation P (prereg §4a) is reported alongside. BH-FDR
across all tested loci, alpha = 0.05. Calibration gate: background pass fraction
must be ~5%.

## Calibration — the null IS calibrated

`[observation]` Fraction of **matched control (background) loci** with matched
empirical P <= 0.05 = **0.044** (263/6,000), i.e. ~5%. The null is therefore
calibrated and loci may be interpreted under the pre-registered rule.

`[observation]` The raw fixed-scale **permutation** P is mildly anti-conservative:
background pass fraction at 0.05 = **0.082** (vs the nominal 0.05). This is the
residual heteroscedasticity of unequal group sizes (G0 n=10,058 vs G1-G5
n=130-679), which an exchangeable permutation of the six group means cannot
absorb. This is why the pre-registered calibrated (matched-control) null, not the
permutation P, is the primary. `[inference]` The permutation P is usable only as a
descriptive cross-check.

## Observed counts and quadrant table

`[observation]` Using `|X| >= 3.16` and matched empirical P <= 0.05:

| quadrant | definition | count | row fraction |
|---|---|---:|---:|
| Q1 | `|X| >= 3.16` and P <= 0.05 | 174 | 0.087 |
| Q2 | `|X| >= 3.16` and P > 0.05 | 1,826 | 0.913 |
| Q3 | `|X| < 3.16` and P <= 0.05 | 263 | 0.044 |
| Q4 | `|X| < 3.16` and P > 0.05 | 5,737 | 0.956 |

`[observation]` **Fraction of low-`|X|` loci passing = 0.044**, identical to the
matched-control null rate (0.044). The high-`|X|` (Akbari-significant) pass
fraction is 0.087 (~2x the null).

`[observation]` **BH-FDR across all 8,000 tested loci: 0 loci pass at FDR <= 0.05
(or <= 0.10).** The smallest matched empirical P is 0.0065.

## Trajectory-level beta for the top loci

`[observation]` Top loci by matched empirical P (all are Akbari-significant,
`|X| >= 3.16`; `beta` is frequency change per year toward the present,
`relatedness_frac = -0.002` for every locus — history time is essentially
`R`-orthogonal here, so relatedness is not inflating these):

| locus | Akbari X | beta | SE | z | concordant steps (of 5) | matched P | BH-FDR |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2:65480598:A:T | -3.27 | -1.45e-5 | 2.1e-6 | -6.93 | 1 | 0.0065 | 0.82 |
| 12:49983378:C:T | -3.16 | -1.44e-5 | 2.3e-6 | -6.20 | 1 | 0.0065 | 0.82 |
| 11:71973305:T:G | -3.79 | -2.43e-5 | 4.0e-6 | -6.03 | 0 | 0.0067 | 0.82 |
| 8:8756191:C:T | +3.31 | +1.63e-5 | 3.0e-6 | +5.47 | 3 | 0.0076 | 0.82 |
| 1:33358762:TGC:T | +3.46 | +2.60e-5 | 5.5e-6 | +4.77 | 4 | 0.0101 | 0.82 |

`[observation]` The strongest-ordering loci (4-5 concordant steps) are all in the
Akbari-significant class, not the low-`|X|` class the hypothesis targets. There is
no low-`|X|` locus with both a strong trend and a concordant ordering that survives
multiple testing.

## Power honesty — stated, not hedged

`[observation]` With m = 6 histories the permutation null has granularity
`1/m! = 1/720 ~ 0.0014`, and the histories are partly nested by descent
(`[observation]` N_eff ~ 4-5; G0 is a documented over-merge). `[inference]` Before
multiple testing, and because 720 permutations cannot resolve a tail beyond
~0.0014, this test **cannot reach genome-wide significance on a 538k-variant
panel by construction**. `[observation]` It did not: **0 loci pass BH-FDR**. The
nominal hits in Q1/Q3 are therefore reported as nominal only and are not
interpreted as genome-wide findings.

## Verdict

`[inference]` **NEGATIVE.** The between-history ordered-rise class is **not
supported** as a locus class that the Akbari correction suppresses:

1. The matched-control null is **calibrated** (background pass fraction 0.044 ~ 5%),
   so the test is valid.
2. Low-`|X|` loci pass at exactly the null rate (0.044 vs 0.044) — **no enrichment**
   of the Akbari-suppressed class. The pre-registered prediction (low-`|X|`
   enrichment above the control ~5%) is **falsified**.
3. The only loci with a strong, correctly-oriented between-history trend are
   Akbari-**significant** (`|X| >= 3.16`, pass fraction 0.087), i.e. loci Akbari
   already flags — the opposite of the hypothesis.
4. **Zero loci** survive BH-FDR, so no locus is even a nominal survivor after
   multiple testing correction.

`[inference]` Combined with the earlier within-trajectory and pooled-replication
lanes, this closes the third, previously untested regime of the over-correction
hypothesis: the between-history ordered rise does not reveal a low-Akbari /
high-replication class under a calibrated, relatedness-aware test.

## Honest limitations

* `[observation]` `relatedness_frac = -0.002`: with these merged groups, history
  time is nearly `R`-orthogonal, so this lane cannot probe the time-GRM collinearity
  regime directly; a finer 13-way grouping (which does show r up to 0.86) would be
  the place to look, but it over-counts replication.
* `[observation]` m = 6 units, N_eff ~ 4-5; the G0 bloc is an over-merge
  (`docs/traj-v3-merge.md`). Effective power is low and the test is
  under-powered for anything but a very strong, very concordant trend.
* `[observation]` The permutation P is anti-conservative (0.082 at nominal 0.05);
  only the matched-control null is trustworthy at locus level.
