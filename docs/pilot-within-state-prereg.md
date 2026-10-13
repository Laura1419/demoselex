# Preregistration: within-state temporal pilot

Registered 2026-10-13, before any within-state computation on genome-wide data.
Companion data: the genome-wide local-ancestry map (`data/derived/elai_genome/`,
`docs/CURRENT-STATE.md` §4). No result-dependent threshold anywhere below.

## Aim

Re-analyse Akbari et al. 2026 genotype data with **local ancestral state as the
demographic correction**, replacing/complementing the GLMM, and check that the two
analyses link up: known selection must survive state-conditioning, carriage must not
be mistaken for selection, and AMY1's position becomes interpretable.

## Hypotheses (fixed now)

- **H1 (sensitivity / "don't lose Akbari's findings")**: variants with Akbari
  `|X| > 5.45` in known-selection loci show |combined within-state z| well above the
  null calibration. If LCT, FADS, SLC45A2 or SLC24A5 go null within-state, the
  conditioning is over-correcting and the project stops for diagnosis.
- **H2 (calibration)**: in null windows, the combined z distribution is centred near 0
  with usable spread; the empirical 99.9th percentile from nulls becomes the full-run
  significance threshold (fixed before the full run).
- **H3 (decomposition, directional expectation made explicit)**: variants whose Akbari
  X is large but within-state z is small are candidate **carriage**; variants large in
  both are candidate **within-background change**. This is the separation the project
  exists to measure. We do NOT predict which variants fall where except for the
  known-selection loci (H1) and AMY1 (measured, no prediction).

## Data

- States: per-window unordered-diploid posterior `q` (1,000 × 55) from
  `data/derived/elai_genome/windows/<tag>/w_s2_k10_seed20261020.q.npy` (seed A;
  seed B used only for the window-level ARI gate).
- Window gate: `ari > 0.5` in `qc_table.csv` (1,922 windows expected). Ungated windows
  are reported but excluded from calibration.
- Individuals: the same 1,000 dated ancients (`samples.json`, `BP` field;
  calendar year = 1950 − BP).
- Genotypes: Akbari genome-wide imputed pgen, hard calls 0/1/2, missing excluded.

## Pilot window set (fixed now, seed 20261020)

- Known-selection: the window containing LCT (chr2:135.85 Mb), FADS (chr11:61.6 Mb),
  SLC45A2 (chr5:33.95 Mb), SLC24A5 (chr15:48.1 Mb) — 4 windows.
- AMY region: chr1 bins 103, 104, 105 (flanks + core-adjacent) — 3 windows.
- Null: 24 windows drawn uniformly from qc_table rows with `ari>0.5` and
  `frac90>0.3`, excluding the above and chr1:100–108 Mb, numpy default_rng(20261020).
- Variants per window: all QC-passing variants in the window, uniform-spaced to ≤600.

## Statistic (per variant, fixed now)

1. Hard-assign each individual to argmax-q state; keep individuals with max q > 0.9.
2. Per state s with n ≥ 60 kept individuals: split dates into 4 quartile bins; require
   ≥ 10 individuals and ≥ 4 called genotypes in every bin, else state dropped
   (variant-level `n_states_used` recorded).
3. Per state: bin alt-allele frequency p_b = sum(g)/(2n_b), binomial SE
   sqrt(p(1−p)/2n_b); WLS slope of p on bin-mean calendar year (per kyr), weights
   1/SE²; z_s = slope/SE.
4. Combined: inverse-variance weighted mean slope across states,
   z_comb = slope_comb / (1/sqrt(sum w_s)). Concordance = weighted fraction of states
   with the dominant sign (weights 1/SE_s²).
5. Variants failing the n-floor are reported as NOT TESTED (not failed).
6. Interpretation guard: within-state slope can legitimately be smaller than the
   marginal slope; that is decomposition, not loss (H3). H1 is judged on the
   known-selection loci only.

## Pilot outputs

`data/derived/pilot_within_state/variants.csv` (one row per variant: window, chrom,
pos, Akbari X if joinable, n_states_used, per-state z/slopes, z_comb, concordance),
`summary.json` (calibration quantiles of null z, known-locus table, AMY table,
X-vs-z correlations), no figures.

## Decision rule

Proceed to the full genome-wide within-state run iff H1 holds for ≥3 of 4
known-selection loci at |z| above the null 99th percentile AND H2 holds (null
|z| median < 1, 99.9th pct finite and recorded). Otherwise diagnose before any
genome-wide run. No thresholds may be revised after seeing results; amendments are
appended and dated, never edited.

---

# AMENDMENT 1 (2026-10-13, before any biological interpretation)

Registered after the first pilot execution (which is reported as a **statistic-validity
failure**, see below) and after discussion of the underlying biology. The hypotheses H1-H3
and the pilot window set stand as originally registered. The *statistic* section is replaced.

## Why the original statistic is void

1. **The binomial z is uncalibrated.** First-run z values reached 1e5-1e6 in NULL windows
   (null 99.9th percentile ~3.2e5 where a calibrated statistic gives ~3). Binomial SEs
   assume independent genomes; individuals sharing a state share history. The calibration
   arm caught this; nothing built on those z values.
2. **The linear slope assumes gradualism.** The registered biological model is punctuated
   adaptation: sharp frequency pulses on specific backgrounds, whose global shadows look
   gradual. A step statistic, not a slope, is the registered shape.
3. **Carriage is not evidence against selection.** Selection may act while a background is
   rare (invisible), surfacing later as apparent carriage. And within-background change
   mixes drift with selection. The decomposition is therefore three-way:
   demography / drift / selection, with drift given an EMPIRICAL baseline, not a formula.

## Revised statistic (replaces the entire Statistic section)

Primary observable: **state enrichment over time**, per (variant, state):

1. **Time bins**: four bins defined ONCE as quartiles of the 1,000 individuals' calendar
   years (1950 - BP); shared across states so enrichment contrasts are simultaneous.
2. **State assignment**: argmax-q at max q > 0.9 (seed A), unchanged.
3. **Enrichment trajectory**: per state s with >= 60 confident individuals, per bin b with
   >= 10 called genotypes: LOR_s(b) = logit(p_s(b)) - logit(p_other(b)), where p_other is
   the allele frequency among all *other confident individuals in the same bin* (all states
   at the same locus, same time). LOR is the enrichment of the allele on background s.
4. **Step statistic** (punctuated model): STEP_s = max over adjacent bin pairs of
   |LOR_s(b+1) - LOR_s(b)|. Also recorded: E0_s = LOR_s at the earliest qualifying bin
   (**enrichment at first detection** - the rare-background-selection signature), and the
   full LOR trajectory.
5. **Background expansion curve**: f_s(b) = fraction of confident individuals in state s
   per bin (the state's own trajectory through time). Timing of enrichment is read
   RELATIVE to f_s: enrichment at/near first detection with flat-or-expanding f_s is the
   registered signature of early selection surfaced by carriage.
6. **Permutation null (within-window, per state)**: permute calendar years among confident
   individuals (preserves state sizes, composition, allele frequencies, date marginals;
   destroys only within-state time structure); recompute STEP_s; 100 permutations in the
   pilot (1,000 in the full run). Per-(variant,state) permutation p-value reported.
7. **Drift baseline (empirical, per state)**: the distribution of STEP across ALL
   (variant, state) pairs in the NULL windows is the drift+noise baseline for this panel;
   thresholds are its 95th/99th/99.9th percentiles, fixed before reading known/AMY arms.
8. **State-specificity (the selection discriminator)**: a (variant, state) pair is a
   selection candidate iff STEP_s exceeds the null baseline 99th percentile AND the max
   STEP among other states at the same variant is below the baseline 50th percentile
   (change on one background, not all). Both numbers recorded per variant.
9. **H1 prediction sharpened**: known-selection loci show enrichment with STEP above the
   drift baseline and permutation p < 0.01, with enrichment present at or before the
   background's expansion (E0 elevated). If a known locus is null under this, diagnose.
10. **Missing data**: genotypes outside 0/1/2 excluded per variant; bins failing the
    n-floor drop that state (NOT TESTED, never failed).

## Decision rule (unchanged in force, restated for the revised statistic)

Proceed to the genome-wide run iff H1 holds for >= 3 of 4 known-selection loci (STEP above
drift baseline q99, permutation p < 0.01) AND H2 holds (null-window STEP distribution
finite, median stable under permutation). No thresholds revised after seeing results.

---

# PILOT RESULT (2026-10-14) — reported per the amended statistic

Execution: 26/31 windows produced data (2 known-selection windows unresolvable: SLC24A5
window c15_w48 has ZERO confident individuals at q>0.9; SLC45A2 window c5_w33 only ~7%;
3 null windows c17_w54/c17_w80/c1_w103-flank also failed their n-floors or produced empty
runs and were excluded). AMY flank windows c1_w104/c1_w105 produced data; c1_w103 did not.
PERMS=100, 30-way then 8-way parallel; two crashed launches from pgen mmap contention
(26 GB file) before the 8-way run completed — recorded as execution notes, not results.

## Drift baseline (24 null windows, 57,600 variant-state pairs)

STEP (max adjacent-bin |Δlogit enrichment|): median 0.0 (most alleles absent/fixed),
q95 = 3.35, q99 = 4.66, q999 = 5.87. Permutation null behaves: 3.1% of null pairs at
perm p < 0.05 (expectation 5% at PERMS=100 resolution with the +1 correction).

## Known-selection loci (H1) — verdict: NOT ESTABLISHED, boundary evidence

- Only 2 of 4 loci produced any confident individuals (LCT c2_w135, FADS c11_w61).
- Their top steps (6.5 at FADS state 19; 6.1 at LCT state 40) sit between q99 (4.66) and
  q999 (5.87): above the 1% tail, below the 0.1% tail.
- But arm-level enrichment is weak: 1.31% of known pairs above q99 vs 1.00% in nulls —
  i.e., the known windows are barely more step-rich than null windows. And the variants
  with the HIGHEST Akbari X in the FADS window (X = 9.1, 8.9) show within-state steps of
  only 0.3-1.2 — the global signal does NOT localize to a single-background step at these
  variants.
- Inspection of the top steps shows the characteristic failure mode of small-n tails:
  p_traj like [0.0, 0.0, 0.093, 0.0] on n=70 (a single carrier appearing in one bin
  produces a LOR swing of ~6). At PERMS=100 the permutation p bottoms at 0.0099. These
  are indistinguishable from carrier-noise at this n.
- Known-window rows carrying Akbari X are sparse (3% join rate — our pilot variant
  sampling spaced 600/window uniformly and missed the tagged variants; a join bug to fix
  in the full run, not a scientific result).

## AMY (measured, no prediction)

c1_w104/c1_w105: 0.89% of pairs above q99 (vs 1.00% null) — no excess of within-state
steps at the AMY flanks beyond background. Top AMY-region steps are again single-carrier
patterns (p_traj [0.015, 0.111, 0.0, 0.143] on n≈84).

## Registered interpretation

1. The statistic is now CALIBRATED (nulls behave), fixing the void first run.
2. H1 is NOT met at face value: with only 2 of 4 loci testable and their top steps below
   the 0.1% baseline tail, we cannot claim known selection survives state-conditioning.
   Per the decision rule this BLOCKS the genome-wide run.
3. The blocker is power, not direction: confident-n per state is 70-130 individuals
   (0-8% of 1000), date bins of ~15-25, and single carriers produce steps as large as
   the top signals. The n-floor (>=60) admitted states whose per-bin n (~15) cannot
   support a step test.
4. Path forward (registered before any further run): (a) join Akbari X for ALL pilot
   variants (the 3% join rate is a sampling artifact); (b) raise the per-bin n-floor so
   a step requires >= 30 called genotypes per bin; (c) test the statistic on an
   INTERNAL positive control with known large n - e.g., pool the two seeds' states or
   use the largest state in each window and confirm the null tail matches; (d) only then
   re-evaluate H1. No thresholds changed retroactively; these are prospective fixes to
   test DESIGN, and the H1 verdict stands as NOT ESTABLISHED at this power.

---

# AMENDMENT 2 (2026-10-14, before the v2 pilot run) — state-space reduction

Motivated by the registered AMENDMENT-1 result: the blocker is per-state n (top states
held 70-130 individuals; per-bin n ~15-25; single carriers produce baseline-sized steps).
Reference-panel enlargement is unavailable. The registered fix is state-space reduction.

## Aggregation rule (fixed now)

1. Per window, from the seed-A q matrix: hard-assign all 1,000 individuals
   (argmax q, no confidence cutoff); rank unordered states by hard count.
2. Keep the **top 5 states** by count; all remaining individuals are pooled as state
   "OTHER". The analysis state space per window is therefore 6 groups.
3. **Soft assignment for counts**: in frequency computations each individual contributes
   weight q_i(state) (its unordered-posterior mass on that state, normalised over the 6
   groups); genotype contribution weighted by q_i * called_i. No q>0.9 cutoff anywhere.
4. Per (variant, group) cell n = sum of weights; a bin enters the test iff
   **weighted n >= 30 called genotypes** (raised from the equivalent of ~15).
5. The four date bins remain the global quartiles of 1950-BP over the 1,000 individuals.
6. seed-B q matrices define the same rule independently; cross-seed reproducibility of
   the aggregated step is reported per window (ARI on the 6 groups).

## Checks before reading the arms

- Composition: NMI(country | 6-group assignment) per window must be >= the 55-state
  median (0.18) for >= 80% of windows, else the merge destroys ancestry signal.
- Null-window drift baseline recomputed under the 6-group rule; the q99 threshold is
  taken from THAT distribution (expected lower than 4.66 because single carriers now
  move group frequencies by ~1/137 rather than ~1/70 - the artifact shrinks with n).

## Expected power

Median top-5 hard mass 0.549 of 1,000 (p10 0.36) -> per-group per-bin weighted n
~31-137, i.e. the floor is met in the median window and in ~90% of windows at bin 1.
A single carrier now shifts a group frequency by <= 1/31 per bin (vs ~1/15 before),
halving the artifact scale, while true pulses retain their full LOR contrast because
enrichment is measured against the OTHER groups pooled.

## Decision rule (unchanged)

Genome-wide run proceeds iff H1 (>= 3 of 4 known-selection loci with step above the
recomputed null q99 and perm p < 0.01; loci failing the confident-n floor count as NOT
TESTED) AND H2 (null calibration holds). Everything else as AMENDMENT 1.

---

# PILOT v2 RESULT (2026-10-14) — AMENDMENT 2 statistic, 31 windows, 200 perms

Execution clean (no failures). 96,000 (variant, group, window) rows; 72,000 null pairs.

## Drift baseline (6-group soft assignment)

q95 = 2.92, q99 = 4.18, q999 = 5.32 (was 3.35/4.66/5.87 under hard q>0.9, 55 states).
Calibration: null frac > q99 = 1.00% (by construction), null frac perm_p<0.05 = 2.68%
(slightly conservative at 200 perms). Single-carrier artifact scale halved as predicted.

## Known-selection loci — verdict: H1 NOT MET (now with power)

- Arm-level: known frac > q99 = 0.64% vs null 1.00% — known windows are NOT more
  step-rich than null windows. Mann-Whitney known>null: p = 1.0; two-sided p < 0.001
  (known upper tail is LOWER: q99 3.95 vs 4.18, q999 5.01 vs 5.32).
- The FADS top step (6.27, group 2|2, perm p 0.005) is now a multi-carrier bump
  (~10 carriers appearing in bin 3 on a ~109-individual group) — but it appears in ONE
  bin and vanishes in the next (p_traj 0 / 0.093 / 0), a spike, not a pulse-with-fixation.
- The LCT group 1|4 trajectory oscillates (0 / 3.2% / 0 / 3.6%) across 4 coarse bins —
  binning resolution, not a temporal signal.
- The canonical LCT variant (rs4988235, chr2:135,851,076) is NOT in our spaced 200-variant
  sample (X-join coverage 1.4% — the known-window variant selection missed tagged
  variants entirely). The 25 LCT-window hits sit at 135.0-136.0 Mb flanks.

## AMY arm

frac > q99 = 0.46% (BELOW null) — no excess of within-state enrichment steps at the AMY
flanks. Top AMY steps (5.7-5.8, groups 2|7, 4|8) again single-bin spikes with perm_p
0.005-0.03, indistinguishable from the null tail.

## Cross-seed reproducibility of the aggregated groups

Median 6-group ARI (seed A vs seed B, same window) = 0.53; mean 0.49; min -0.05.
Aggregation is reproducible in most windows but not all.

## Registered interpretation

1. H1 is NOT MET, and at this power the verdict is informative: the within-state step
   statistic, calibrated on 72k null pairs, does NOT show the known-selection loci as
   outliers above background. The direction of the deficit (known < null upper tail)
   suggests the coarse 4-bin grouping AVERAGES OVER the pulses the model predicts.
2. Two design defects remain before any genome-wide run: (a) variant selection must
   include ALL Akbari-X-tagged variants per window (join coverage was 1.4%); (b) 4
   genome-wide quartile bins spanning -11,483..1800 AD compress Neolithic-transition
   pulses (the interesting interval) into one or two bins. Bins must be era-adaptive
   (e.g., fixed boundaries at -8500, -5000, -3000, -1500, 0, 1800 AD chosen BEFORE
   looking at any allele trajectory).
3. Registered next step: pilot v3 with (a) + (b) on the same 31 windows. H1 remains
   NOT ESTABLISHED; no genome-wide run.
