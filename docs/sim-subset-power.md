# Simulation: the subset-selection regime and the max-over-subsets device

Produced by `sim/sweeps/run_subset_power.py` on `sim/scenarios/subset_selection.yaml`.
Raw table: `sim/sweeps/results/subset_power.csv` (42 rows). Figure: `fig/subset_power.png`.

Design: `K = 6` trajectories (matching the order of magnitude of the independent
histories actually inferred from AADR data). A fraction `f` of trajectories experience
selection on the tested locus (`f = 0.2, 0.4, 0.6, 1.0`, so 1, 2, 4 and 6 of 6), with
selection strength `0.5, 1.0, 2.0`. Three statistics are compared, all calibrated
against the neutral baseline at the 95th percentile:

- `akbari` - the reconstructed structure-corrected GLMM;
- `traj_all` - concordance/meta-analysis over all trajectories;
- `traj_max_subset` - the **max-over-subsets** statistic (maximum |z| over subsets of
  size 2-4), the pre-registered statistic in `docs/subset-regime-prereg.md`.

## Calibrated power

| f (selected/K) | strength | akbari | traj_all | **max_subset** |
|---|---:|---:|---:|---:|
| 0.2 (1/6) | 0.5 | 0.20 | 0.20 | **0.45** |
| 0.2 (1/6) | 1.0 | 0.25 | 0.30 | **0.95** |
| 0.2 (1/6) | 2.0 | 0.30 | 0.30 | **1.00** |
| 0.4 (2/6) | 0.5 | 0.60 | 0.50 | **0.85** |
| 0.4 (2/6) | 1.0 | 0.80 | 0.75 | **1.00** |
| 0.6 (4/6) | 0.5 | 1.00 | 0.90 | 1.00 |
| 1.0 (6/6) | any | 1.00 | 1.00 | 1.00 |

## Calibrated FPR

| scenario | akbari | traj_all | max_subset |
|---|---:|---:|---:|
| neutral_drift | 0.050 | 0.050 | 0.050 |
| replacement_only | 0.045 | 0.071 | **0.112** |
| subset_selection cells | 0.03-0.08 | 0.03-0.06 | 0.04-0.06 |

## The uncalibrated device is worthless - this is the key caveat

`fpr_naive_device` for `traj_max_subset` is **0.55-0.64** in every cell, versus
0.04-0.06 once the max is calibrated against the *same* max-over-subsets statistic on
matched controls. A max over ~50 subsets is a ~50-fold multiple-testing device: comparing
it to a non-max null would fabricate discoveries at a ~60% rate. Any real-data use of
this statistic **must** carry its own max-over-subsets null.

## Answer to the central question

**[observation]** Where selection acts in a **minority** of trajectories (`f ~ 0.2`),
the max-over-subsets statistic has calibrated power **0.45-1.00** against the GLMM's
**0.20-0.30** - a 2-4x gain, and the gain is largest at moderate effect size. This is
the one regime in which a trajectory method beats the pooled GLMM, and it is now
quantified rather than asserted.

**[observation]** As `f` grows the advantage disappears: by `f = 0.6` every statistic
has power ~1.0, and the GLMM is as good as anything.

## Limits and the price

- **FPR cost.** Under pure single-source replacement, the calibrated max-over-subsets
  FPR is **0.112** versus Akbari's **0.045** - the same anti-conservative failure mode
  documented in `docs/G-validation-report.md` (Decision 3). The device is therefore
  acceptable only where the matched null absorbs replacement, which is exactly the
  assumption the Tier-2 analysis must stress.
- **Argmax behaviour.** The best subset is the maximum allowed size (4) for ~87% of
  loci, including controls; the device mostly inflates the statistic rather than
  identifying a genuinely small agreeing subset. A pre-specified, smaller subset or a
  complement-flatness score would be better behaved and is a recommended follow-up.
- 6 trajectories, simple drift model, no linkage, no trajectory-inference error beyond
  the earlier sweeps. Absolute power values are optimistic for real data.

## Implication for the real-data test

This simulation says precisely where to look: loci where a **minority** of the
independent histories move together while the rest are flat. The pre-registered
real-data test of that regime (`docs/subset-regime-results.md`) returned
**0 of 2,000 candidates at BH-FDR 0.05** (minimum q = 0.639) with candidate enrichment
of only 1.13x over the matched max-null. So at Tier-1 power the regime is not
detectable - which is consistent with the simulation's requirement of effect sizes at
or above `strength 1.0` to reach power ~0.95 with a clean 1-of-6 pattern.
