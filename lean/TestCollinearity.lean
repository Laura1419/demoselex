import TrajectorySelection.Collinearity

namespace TrajectorySelection.Collinearity

-- A unit GRM mode along the first of two observations.
private def mode : Fin 2 → ℝ := ![1, 0]
private def time : Fin 2 → ℝ := ![3, 4]
private def alignedTime : Fin 2 → ℝ := ![3, 0]
private def nearTime : Fin 2 → ℝ := ![3, 1]

example : energy mode = 1 := by norm_num [energy, dot, mode, Fin.sum_univ_two]
example : energy time = 25 := by norm_num [energy, dot, time, Fin.sum_univ_two]
example : cosineSquared time mode = 9 / 25 := by
  norm_num [cosineSquared, energy, dot, time, mode, Fin.sum_univ_two]
example : energy (residual time mode) = 16 := by
  norm_num [energy, residual, dot, time, mode, Fin.sum_univ_two]
example : projectedVIF time mode = 25 / 16 := by
  norm_num [projectedVIF, energy, residual, dot, time, mode, Fin.sum_univ_two]
-- Positive finite covariance retains information even for perfectly aligned time.
example : workingInformation alignedTime mode 1 8 = 1 := by
  norm_num [workingInformation, energy, residual, dot, alignedTime, mode,
    Fin.sum_univ_two]

-- At 90% cosine alignment and a mode/noise ratio of 99, at most
-- 10% + 1% of the isotropic time information remains.
example : workingInformation nearTime mode 1 99 / (energy nearTime / 1) ≤
    (1 / 10 : ℝ) + 1 / 100 := by
  apply near_alignment_working_bound nearTime mode 1 99 (1 / 10) (1 / 100)
  all_goals norm_num [energy, dot, nearTime, mode, cosineSquared, Fin.sum_univ_two]

-- Inspect inherited trust as well as the absence of project-level sorry/axioms.
#print axioms TrajectorySelection.Collinearity.residual_energy_cosine
#print axioms TrajectorySelection.Collinearity.projectedVIF_eq
#print axioms TrajectorySelection.Collinearity.working_information_fraction
#print axioms TrajectorySelection.Collinearity.near_alignment_working_bound

end TrajectorySelection.Collinearity
