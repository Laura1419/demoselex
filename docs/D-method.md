# D. Proposed method

> Owner task: lean-* + sim-stat. Status: TODO.

## Definitions

The formal definitions are in `lean/TrajectorySelection/Defs.lean` (imported by
`lean/TrajectorySelection.lean`). For an index set `I` and population-state type `P`,
a temporal trajectory consists of `population : I → P` and `time : I → ℝ`.
An allele-frequency path `p : I → ℝ` has `0 ≤ pᵢ ≤ 1` at every index. The
index can represent dated samples or aggregate time bins. These definitions
alone impose neither ordered/unique times nor ancestry independence; those
are assumptions to state in later models, not consequences of the data type.

For a finite observed set `S` and integer sample sizes `nᵢ`, let
`n = Σᵢ∈S nᵢ`, `t̄ = Σ nᵢtᵢ / n` and `p̄ = Σ nᵢpᵢ / n`. The sample-size-weighted
least-squares slope is

`β̂ = Σ nᵢ(tᵢ − t̄)(pᵢ − p̄) / Σ nᵢ(tᵢ − t̄)²`.

Lean returns `none` if `n = 0` or the weighted time variation is zero; it
never silently treats an undefined slope as zero. `directionSign` maps a
slope to `−1`, `0`, or `1` by exact sign. For a finite set of trajectory slopes,
let `N₊` and `N₋` count strictly positive and negative slopes; directional
concordance is `max(N₊, N₋)/(N₊ + N₋)`, or `none` if both counts vanish.
Zero slopes are excluded, and a tie gives `1/2`. This descriptive statistic
is **not** a calibrated p-value, evidence of selection, or a claim of
independent histories.

In the replacement-only model, two source frequencies `pA,pB ∈ [0,1]`
are constant while mixing weights `wᵢ ∈ [0,1]` vary. The pooled path is
`pᵢ = (1 − wᵢ)pA + wᵢpB`. Thus a pooled slope can arise from turnover
alone. `TrajectoryCovariance` records a symmetric, positive-semidefinite
finite-support covariance kernel on trajectory labels. Its off-diagonal
entries can represent shared history; no normalization, effective replicate
count, likelihood, or covariance estimator is asserted here. Downstream
formal results must add assumptions for null exchangeability, independence,
noise and identifiability explicitly.

The inherited `Core.lean`, `FrequencyPath.lean` and `Concordance.lean`
scaffolds have older convenience definitions and provisional statistical
placeholders; see `lean/AXIOMS.md`. In particular their `check_replication`
is **not** the directional-concordance definition above. No biological
inference should use those placeholders as calibrated statistics.

## Stage A: trajectory inference
## Stage B: replicated directional change
## Statistics (Model 0 .. Model 3)
## Null calibration
## One-paragraph plain-language summary
> "We inferred major population trajectories from genome-wide variation independently
> of each tested locus. We then estimated allele-frequency change through time along
> each trajectory and tested whether the direction and magnitude of change were
> concordant across genetically distinct trajectories, while accounting for their
> non-independence. Significance was calibrated against matched genomic control loci
> and demographic simulations."
