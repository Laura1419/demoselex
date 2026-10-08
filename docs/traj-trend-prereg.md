# Pre-registration — between-history ordered rise (traj-trend)

> Owner lane: traj-trend. **Written BEFORE any result was produced.** The
> timestamped ordering is enforced by this task's procedure: `analysis/traj_trend.py`
> is implemented and run only *after* this file exists, and
> `docs/traj-trend-results.md` is written only *after* the run. No result,
> count, or P-value in this file is anything but a pre-specified expectation.

## 1. Hypothesis (third regime, distinct from the two tested before)

The two earlier lanes tested *within-trajectory* change (a trajectory rising over
its own time span) and pooled replication across trajectories. This lane tests a
different, previously untested pattern:

> **Between-history ordered rise.** Every individual history is internally *flat*,
> but each successive history — ordered in time — carries a **higher** mean allele
> frequency.

`[speculation]` If selection acts on a standing allele along a shared,
time-ordered ancestry that is later sampled as several distinct histories, the
locus can rise monotonically across histories while being flat *within* each one.
A within-trajectory slope sees nothing; the pooled `MVN(0, sigma^2 K)` random
effect in the Akbari GLMM absorbs a time-tracking between-group difference as
structure. This is therefore the best remaining candidate for a locus class the
Akbari correction could suppress.

## 2. Unit of analysis (decided before the run)

**Unit = `merged_group` (G0..G5), m up to 6 histories.**

Justification, from `docs/traj-v3-merge.md`:
* The 13 `trajectory_id` values are known to be over-split by sampling period and
  admixture fraction (e.g. T1/T3/T6 Avar split, T0/T4 medieval split, pairwise r up
  to 0.86). Treating them as 13 independent histories double-counts shared descent.
* `merged_group` is the repository's defensible replication unit; it is the label
  that `data/derived/traj_relatedness_v3.tsv` is indexed by, so unit and relatedness
  matrix match with no hand-set values.
* `group_id` is the AADR *published* annotation and is explicitly **not** used.

Honest caveats carried into the run:
* `merged_group` is itself an over-merge at G0 (a broad west/central-Eurasian bloc,
  `[observation]` in the merge doc). The effective count of genuinely distinct
  histories is N_eff ~ 4-5, not 6.
* `unassigned` individuals (`trajectory_id == -1`) are dropped.

Relatedness `R` = the 6x6 `data/derived/traj_relatedness_v3.tsv` matrix, used
directly (it already aggregates the 13-way matrix to merged groups and is verified
positive-definite by `sim.trajectory_stat._correlation`). No further aggregation.

## 3. Statistic (exact formula)

For each locus, let `p_i` be the mean ALT allele frequency of history `i`
(`p_i = mean over individuals in group i of dosage/2`) and let `t_i` be its
temporal placement, `t_i = - mean(date_BP)` over the group so that **larger `t` =
more recent**. Fit, at the level of the m histories:

```
p = alpha * 1 + beta * t + e ,   e ~ N(0, sigma^2 R)
```

**GLS** (center by the GLS-weighted means):
```
w        = R^{-1} 1
a        = (1' R^{-1} t) / (1' R^{-1} 1)
t_c      = t - a
p_c      = p - (1' R^{-1} p) / (1' R^{-1} 1)
beta_hat = (t_c' R^{-1} t_c)^{-1} (t_c' R^{-1} p_c)
sigma^2  = (r' R^{-1} r) / (m - 2) ,   r = p_c - beta_hat * t_c
SE(beta) = sqrt( sigma^2 / (t_c' R^{-1} t_c) )
z        = beta_hat / SE(beta)
```

Additional reported quantities (same definition for every locus and every control):
* `m` — number of histories present for the locus.
* `ordering_concordant` — number of successive increases in `p` when histories are
  sorted by ascending `t` (range 0..m-1). This is the pattern the hypothesis predicts.
* `relatedness_frac` — fraction of temporal leverage absorbed by relatedness,
  `1 - (t_c' R^{-1} t_c) / (t_c' t_c)`. `0` means `t` is orthogonal-in-`R` to the
  independent design; larger means relatedness and time are collinear here.

Sign convention: `beta_hat > 0` = frequency **rises** toward the present, the
predicted direction.

## 4. Nulls

**(a) Permutation of `p` across histories.** Holding `R` and `t` fixed, permute the
`p_i` across the m histories (all `m! = 720` permutations for m=6) and recompute
`z`. Permutation P `= (1 + #{|z_perm| >= |z_obs|}) / (1 + n_perm)`. This is the
primary P.

**(b) Matched control loci in the genome.** Background panel loci with
`|X| < 3.16`, matched to a candidate on chromosome, Akbari `AF` (±0.10 caliper,
widening 1x/2x/4x), missingness and temporal coverage, same-chromosome and >= 1 Mb
apart. Matched empirical P `= (1 + #{|z_ctrl| >= |z_cand|}) / (1 + n_ctrl)`.
This is the calibration null.

## 5. Decision rule (fixed before the run)

1. Primary empirical P from null (a) for every panel locus.
2. Benjamini-Hochberg FDR across all tested panel loci; **alpha = 0.05**.
3. Report the quadrant counts at `|X| >= 3.16` vs `P <= 0.05`, the fraction of
   low-`|X|` loci passing, and the matched-control pass fraction.

**Calibration gate:** the fraction of *matched control* (background) loci with
`P <= 0.05` must be ~5%. If it is not, the null is mis-calibrated and **no locus
may be interpreted** until it is fixed.

## 6. Panel cap (fixed before the run)

Use the prepared cache `data/derived/tier2_pool_cache.pkl` (538,274 QC-confirmed
Akbari variants, columns incl. `X`, `AF`, `is_candidate`, `pgen_index`), copied from
`/tmp/tier2_pool_cache.pkl` so it survives `/tmp` cleanup. Default run cap
**8,000 variants** (< 12 min, within the 5,000-20,000 band); the Validation run uses
`--max-variants 2000`. Up to `cap/4` candidates (`|X| >= 3.16`) and the remainder
background (`|X| < 3.16`). Exactly what was analysed is reported in the results doc.

## 7. Pre-registered prediction and falsification

`[speculation]` **Prediction.** If a between-history ordered-rise class exists, a few
low-`|X|` loci will show a positive `beta_hat` with a concordant ordering and a small
permutation P, at a rate above the matched-control ~5%.

**Falsifiers (any of these is a valid negative result):**
* The log of low-`|X|` loci passing equals the matched-control pass fraction
  (no enrichment), or the null is mis-calibrated.
* `relatedness_frac` is high for the top loci — i.e. the apparent "trend" is exactly
  the time-history collinearity the GLMM is designed to absorb, not signal.
* `beta_hat` signs are not concordant across the history ordering.

## 8. Power honesty (stated before the run)

With `m` up to 6 histories the permutation null has granularity `1/m!` = 1/720
~ 0.0014 before multiple testing, and the histories are partly nested by descent
(N_eff ~ 4-5). Genome-wide significance over the 538,274-variant panel requires
`P` far below what 720 permutations can resolve, so **by construction this test
cannot reach genome-wide significance on the full panel**; only nominal/regional
statements are possible, and they are reported as such, not as genome-wide hits.
