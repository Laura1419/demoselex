import TrajectorySelection.Replacement

open TrajectorySelection

-- The statistic uses the existing definition, not a hand-coded event indicator.
example : pairConcordance (true, true) = some 1 := by
  norm_num [pairConcordance, directionalConcordance, signSlope, directionSign,
    show (Finset.univ : Finset (Fin 2)) = {0, 1} by decide,
    Finset.filter_insert, Finset.filter_singleton]

example : pairConcordance (true, false) = some (1 / 2 : ℝ) := by
  norm_num [pairConcordance, directionalConcordance, signSlope, directionSign,
    show (Finset.univ : Finset (Fin 2)) = {0, 1} by decide,
    Finset.filter_insert, Finset.filter_singleton]

-- An independent biased sign model (q = 3/4) has strictly higher mass at 1.
example : highConcordanceMass (3 / 4) = (5 / 8 : ℝ) := by
  rw [highConcordanceMass_eq]
  norm_num

-- Different replacement histories cannot alter concordance for fixed observed signs.
example {Time : Type*} (a b : ReplacementOnly Time) :
    (ReplacementScenario.mk a (true, false)).concordance =
      (ReplacementScenario.mk b (true, false)).concordance :=
  replacement_concordance_determined_by_signs a b (true, false)

-- Without sign bias, the high-concordance mass need not exceed the null.
example : highConcordanceMass (1 / 2) = (1 / 2 : ℝ) := replacement_high_mass

-- Even q = 0 (perfect shared negative signs) gives a high-concordance event:
-- this is a two-sided statistic, not a directional positive-only test.
example : highConcordanceMass 0 = 1 := by
  rw [highConcordanceMass_eq]
  norm_num
