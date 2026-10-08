# Pre-registration: subset-selection regime test (max-over-subsets)

Status: **WRITTEN BEFORE RESULTS.** This file was authored before
`analysis/subset_regime.py` was run and before `data/derived/subset_regime.tsv`
existed. It is git-untracked at the time of writing (new file, no prior commit);
the parent commits it. No result from the subset scan was inspected before this
text was fixed. Any later deviation is recorded explicitly in
`docs/subset-regime-results.md` and the affected analysis is marked exploratory.

## Motivation (one regime only)

The first-pass genome-wide comparison found **no** low-|X|/high-replication
excess (`docs/I-quadrant-first-pass.md`). The reason is structural: Akbari's GLMM
fits **one** time slope across all individuals, so if selection acts in only a
**subset** of the population histories while the rest stay flat or move the other
way, the pooled slope is diluted and / but the per-subset concordance is high. This
is the single regime in which a trajectory/within-group method is expected to beat
the pooled GLMM, and sim-sweeps already hinted at it (`one_population_selection`
power 0.667 vs Akbari 0.167). This pre-registration is a targeted test of exactly
that regime. It is **not** a re-test of the pooled first pass.

## (a) Replication units

`data/derived/trajectories_v3.tsv` **does not exist** (verified before running).
Therefore the replication units are the **13 trajectories** in
`data/derived/trajectories.tsv` (`trajectory_id` 0..12; the `-1` unassigned
individuals are dropped, as in `analysis/reanalysis_scan.py`). We use the merged
v3 grouping if it had existed; it does not, so we use the 13-trajectory v2
grouping. Each trajectory is its own independent group `G`; the prespecified
relatedness matrix `data/derived/traj_relatedness.tsv` is applied inside the GLS
meta-analysis exactly as in `analysis/reanalysis_scan.py` (working independence
would be the fallback only if that file were absent).

## (b) Statistic (max-over-subsets)

For each locus, let `G` be the set of trajectories with the required support
(>= 2 distinct time bins, as in `analysis/reanalysis_scan.py`). Enumerate every
subset `S` of `G` with **2 <= |S| <= 4**. For each `S`, compute the sample-size
weighted frequency slope and binomial SE for each trajectory in `S` (the
estimator in `sim/trajectory_stat.py`, `_weighted_slope`), then the correlated
fixed-effect GLS meta combined slope and its `z = meta_slope / meta_se` restricted
to `S` (same GLS math as `trajectory_statistic`, implemented on the precomputed
per-trajectory slopes/SEs). The **locus statistic** is

    T(locus) = max over admissible S of |z(S)|.

We also record: the argmax subset size `|S*|`, the signed `z(S*)`, the meta slope
of the **complement** `G \ S*` when `|G \ S*| >= 2` (else NaN), and two boolean
flags: `complement_flat` (the complement meta slope is consistent with zero) and
`complement_opposite` (the complement meta slope has the opposite sign of
`z(S*)`). These flags characterise whether the best subset separates a coherent
sub-block from a flat/opposing remainder - the signature of the target regime.

This max-over-subsets device is deliberately a multiple-testing machine; it is
only valid because it is calibrated against the **same** statistic on controls
(see (c)).

## (c) Null: matched control AADR sites, same max-over-subsets statistic

Controls are the AADR sites in the same panel that are **not** Akbari-significant
(`|X| < 3.16`), matched to each candidate with the caliper matcher of
`analysis/reanalysis_scan.py`: same chromosome, `|pos difference| >= 1 Mb`
exclusion, and calipers on Akbari allele frequency, missingness and temporal
coverage, widening 1x/2x/4x until at least `--min-controls` (default 20) survive;
last-resort fallback is all same-chromosome background loci. Background loci are
excluded from the candidate set.

Crucially, `T` is computed for **every** panel locus (candidates and background)
by the identical max-over-subsets procedure. The empirical upper-tail P for a
candidate is `(1 + #{controls with T_control >= T_candidate}) / (m + 1)`. The
control statistic is the **max-over-subsets** statistic, so the calibration
absorbs the selection over subsets.

## (d) Decision rule

* Empirical upper-tail P per candidate (formula above).
* **BH-FDR** across the tested candidates at **alpha = 0.05**.
* **Minimum subset size = 2** (maximum subset size = 4).
* A locus "passes" if its BH-FDR q <= 0.05. `|X| < 3.16` defines "low-|X|";
  `|X| >= 5.45` is the genome-wide Akbari cut; `3.16 <= |X| < 5.45` is the
  intermediate band.
* A candidate with `T >= 3.16`-equivalent support but below the genome-wide cut
  is only interesting if its BH-FDR q <= 0.05 **and** its argmax subset is a
  strict subset (`|S*| < |G|`) with a flat/opposite complement.

## (e) Pre-registered prediction, and what would falsify it

Prediction: in the subset-selection regime, at least some Akbari-significant
(`|X| >= 3.16`) loci - in particular loci with only **modest** |X|
(`3.16 <= |X| < 5.45`) - show a max-over-subsets `T` that (i) exceeds their
matched-control null and (ii) attains BH-FDR q <= 0.05, with the winning subset
being a **proper** subset whose complement is flat or opposite. Equivalently, the
candidate mean `T` exceeds the matched-control mean `T` at these loci.

Falsification: the prediction is **falsified** if (a) no candidate passes
BH-FDR at alpha 0.05, or (b) candidates show no enrichment of `T` over matched
controls (candidate-mean `T` <= control-mean `T`), or (c) every passing candidate
has a full-set argmax `|S*| = |G|` (i.e. the signal is not subset-restricted, so
the GLMM would have seen it too). Any of these is a valid, reportable negative
result and must be reported plainly.

## Honest limitations fixed in advance

* Tier-1 panel: ~26% of Akbari's individuals; low-coverage pseudo-haploid
  genotypes; both statistics attenuated.
* Trajectories are partly over-split / partly stable (see
  `docs/I-quadrant-first-pass.md`), so `|G|` and the relatedness matrix are
  approximations.
* The max-over-subsets statistic is anti-conservative by construction; only the
  matched max-over-subsets null is used, and results are still Tier-1 scouting,
  not a discovery list.

## Implementation note (fixed before running)

`analysis/subset_regime.py` reuses `analysis/reanalysis_scan.py` for cohort
loading, the AADR `.snp`/Akbari join, the byte-random-access genotype reader and
the panel builder (SOURCES), then computes `T` on the precomputed per-trajectory
slopes/SEs (no re-fitting), writes `data/derived/subset_regime.tsv` incrementally
after each genotype chunk, and reports counts against the rule in (d)/(e).
