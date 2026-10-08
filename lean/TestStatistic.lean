import TrajectorySelection.Statistic

open TrajectorySelection

-- Undefined at zero informative trajectories, tied counts, unanimity.
example : concordanceCounts 0 0 = none := by norm_num [concordanceCounts]
example : concordanceCounts 1 1 = some (1 / 2 : ℝ) := by norm_num [concordanceCounts]
example : concordanceCounts 2 0 = some (1 : ℝ) := by norm_num [concordanceCounts]

-- At n=2, both unanimous outcomes have C=1; ties have C=1/2.
example : fairNullConcordance 2 0 = 1 := by norm_num [fairNullConcordance]
example : fairNullConcordance 2 1 = 1 / 2 := by norm_num [fairNullConcordance]
example : fairNullMean 2 = 3 / 4 := by
  norm_num [fairNullMean, fairBinomialAverage, fairBinomialWeight,
    fairNullConcordance, Finset.sum_range_succ]
example : fairNullVariance 2 = 1 / 16 := by
  norm_num [fairNullVariance, fairNullMean, fairBinomialAverage,
    fairBinomialWeight, fairNullConcordance, Finset.sum_range_succ]

-- The new label 2 is a near-copy of label 0: both map to family 0.
example (ρ : ℕ → ℕ → ℝ) :
    effectiveTrajectoryCount (insert 2 ({0, 1} : Finset ℕ))
        (fun i => if i = 2 then 0 else i) ρ =
      effectiveTrajectoryCount ({0, 1} : Finset ℕ)
        (fun i => if i = 2 then 0 else i) ρ := by
  apply effectiveTrajectoryCount_duplicate
  simp

#print axioms concordanceCounts_range
#print axioms effectiveTrajectoryCount_duplicate
#print axioms effectiveFamilyCount_range
#print axioms fairNullVariance_formula
