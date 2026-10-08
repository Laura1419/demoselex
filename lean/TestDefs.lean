import TrajectorySelection.Defs

open TrajectorySelection

-- Exact signs distinguish negative, zero, and positive slopes.
example : directionSign (-1 : ℝ) = -1 := by simp [directionSign]
example : directionSign (0 : ℝ) = 0 := by simp [directionSign]
example : directionSign (1 : ℝ) = 1 := by simp [directionSign]

-- A single positive trajectory has full directional concordance.
example : directionalConcordance ({0} : Finset ℕ) (fun _ => (1 : ℝ)) = some 1 := by
  simp [directionalConcordance, directionSign]

-- One zero slope has no directional concordance denominator.
example : directionalConcordance ({0} : Finset ℕ) (fun _ => (0 : ℝ)) = none := by
  simp [directionalConcordance, directionSign]

-- A positive/negative tie counts half of the nonzero trajectories.
example : directionalConcordance ({0, 1} : Finset ℕ)
    (fun i => if i = 0 then (1 : ℝ) else -1) = some (1 / 2 : ℝ) := by
  norm_num [directionalConcordance, directionSign, Finset.filter_insert, Finset.filter_singleton]

-- Zero total sample size never produces a slope.
example {Index Population : Type*} (trajectory : TemporalTrajectory Index Population)
    (path : AlleleFrequencyPath Index) (observed : Finset Index) :
    sampleWeightedSlope trajectory path observed (fun _ => 0) = none := by
  simp [sampleWeightedSlope, weightedMeanTime, totalSampleSize]

-- Identical source frequencies cannot make the pooled path move by turnover alone.
example {Index : Type*} (m : ReplacementOnly Index) (i : Index)
    (h : m.sourceA = m.sourceB) : m.pooled i = m.sourceB := by
  simp [ReplacementOnly.pooled, h, sub_mul]
