import TrajectorySelection.Collinearity

namespace TrajectorySelection.Collinearity

-- A unit GRM mode along the first of two observations.
private def mode : Fin 2 → ℝ := ![1, 0]
private def time : Fin 2 → ℝ := ![3, 4]
private def alignedTime : Fin 2 → ℝ := ![3, 0]

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

-- Inspect inherited trust as well as the absence of project-level sorry/axioms.
#print axioms TrajectorySelection.Collinearity.residual_energy_cosine
#print axioms TrajectorySelection.Collinearity.projectedVIF_eq
#print axioms TrajectorySelection.Collinearity.working_information_fraction

end TrajectorySelection.Collinearity
