import Mathlib.Basic.Real.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic

namespace TrajectorySelection

/- A history is an indexed family of population states, each observed at a real-valued time. -/
/-- A temporal population trajectory separates the sample index from the population state and its time. -/
structure TemporalTrajectory (Index Population : Type*) where
  population : Index → Population
  time : Index → ℝ

/- A locus has one frequency at each index of a trajectory; the bounds express a proportion. -/
/-- An allele-frequency path on the indices of a temporal trajectory. -/
structure AlleleFrequencyPath (Index : Type*) where
  frequency : Index → ℝ
  bounded : ∀ i, 0 ≤ frequency i ∧ frequency i ≤ 1

/- The total sample size is the weight available for a finite set of observations. -/
/-- Sum of the nonnegative integer sample sizes in a finite set of indices. -/
def totalSampleSize {Index : Type*} (observed : Finset Index) (sampleSize : Index → ℕ) : ℝ :=
  ∑ i ∈ observed, (sampleSize i : ℝ)

/- A weighted average of times exists only if at least one observation has a positive sample size. -/
/-- Sample-size-weighted mean time, undefined for zero total sample size. -/
noncomputable def weightedMeanTime {Index Population : Type*}
    (trajectory : TemporalTrajectory Index Population) (observed : Finset Index)
    (sampleSize : Index → ℕ) : Option ℝ :=
  let weight := totalSampleSize observed sampleSize
  if weight = 0 then none
  else some ((∑ i ∈ observed, (sampleSize i : ℝ) * trajectory.time i) / weight)

/- A weighted average of frequencies uses the same samples and weights as the time mean. -/
/-- Sample-size-weighted mean allele frequency, undefined for zero total sample size. -/
noncomputable def weightedMeanFrequency {Index : Type*}
    (path : AlleleFrequencyPath Index) (observed : Finset Index)
    (sampleSize : Index → ℕ) : Option ℝ :=
  let weight := totalSampleSize observed sampleSize
  if weight = 0 then none
  else some ((∑ i ∈ observed, (sampleSize i : ℝ) * path.frequency i) / weight)

/- This is weighted least-squares frequency change per unit time, not a significance estimate. -/
/-- Weighted OLS slope; returns `none` when there is no weight or no weighted time variation. -/
noncomputable def sampleWeightedSlope {Index Population : Type*}
    (trajectory : TemporalTrajectory Index Population) (path : AlleleFrequencyPath Index)
    (observed : Finset Index) (sampleSize : Index → ℕ) : Option ℝ := do
  let meanTime ← weightedMeanTime trajectory observed sampleSize
  let meanFrequency ← weightedMeanFrequency path observed sampleSize
  let variance := ∑ i ∈ observed,
    (sampleSize i : ℝ) * (trajectory.time i - meanTime) ^ 2
  if variance = 0 then none
  else some ((∑ i ∈ observed,
    (sampleSize i : ℝ) * (trajectory.time i - meanTime) *
      (path.frequency i - meanFrequency)) / variance)

/- The integer sign records rising, falling, or exactly flat change without a cutoff. -/
/-- Sign of a slope: `1` for positive, `-1` for negative, and `0` for zero. -/
noncomputable def directionSign (slope : ℝ) : ℤ :=
  if 0 < slope then 1 else if slope < 0 then -1 else 0

/- Flat slopes are excluded; a tied positive/negative split has concordance one half. -/
/-- Fraction of nonzero slopes in the more common direction; `none` if all slopes are zero. -/
noncomputable def directionalConcordance {Trajectory : Type*}
    (selected : Finset Trajectory) (slope : Trajectory → ℝ) : Option ℝ :=
  let positive := (selected.filter fun t => directionSign (slope t) = 1).card
  let negative := (selected.filter fun t => directionSign (slope t) = -1).card
  if positive + negative = 0 then none
  else some ((max positive negative : ℕ) / (positive + negative : ℕ) : ℝ)

/- Turnover can change a pooled frequency without changing either source's own frequency. -/
/-- Two constant-frequency sources with a time-varying mixing fraction in the unit interval. -/
structure ReplacementOnly (Index : Type*) where
  sourceA : ℝ
  sourceB : ℝ
  sourceA_bounded : 0 ≤ sourceA ∧ sourceA ≤ 1
  sourceB_bounded : 0 ≤ sourceB ∧ sourceB ≤ 1
  mixing : Index → ℝ
  mixing_bounded : ∀ i, 0 ≤ mixing i ∧ mixing i ≤ 1

/- At each index, the pooled allele frequency is the convex mixture of fixed sources. -/
/-- Pooled frequency path under replacement only, with no within-source change. -/
def ReplacementOnly.pooled {Index : Type*} (model : ReplacementOnly Index) (i : Index) : ℝ :=
  (1 - model.mixing i) * model.sourceA + model.mixing i * model.sourceB

/- Shared history may correlate trajectories; a covariance must be symmetric and nonnegative in quadratic form. -/
/-- A finite-support positive-semidefinite covariance kernel over trajectory labels. -/
structure TrajectoryCovariance (Trajectory : Type*) where
  cov : Trajectory → Trajectory → ℝ
  symmetric : ∀ a b, cov a b = cov b a
  positiveSemidefinite : ∀ (selected : Finset Trajectory) (coeff : Trajectory → ℝ),
    0 ≤ ∑ i ∈ selected, ∑ j ∈ selected, coeff i * cov i j * coeff j

end TrajectorySelection
