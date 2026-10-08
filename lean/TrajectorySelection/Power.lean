/-
Power separation between replacement and replicated selection

AXIOMS policy: Every `sorry`, `axiom`, and `native_decide` must be documented in lean/AXIOMS.md
-/

import TrajectorySelection.Defs
import TrajectorySelection.Statistic
import TrajectorySelection.Replacement
import Mathlib.Tactic

namespace TrajectorySelection

/-!
# Power separation

This file formalizes the finite, sign-level comparison between *replacement-only*
selection (an exchangeable-sign null) and *replicated* directional selection. It builds
on the definitions in `Defs.lean` (directional concordance, `ReplacementOnly`),
`Statistic.lean` (the fair folded-binomial null, effective trajectory count) and
`Replacement.lean` (the two-sign scenario).

The probability models below are explicit **finite** sign models. Statements of the form
"directional selection has power" are statements about these declared masses, not about
how a physical drift/noise ratio produces them.
-/

/-- Slope vector of a sign function: `+1` at positive signs and `-1` at negative signs.
    A nonzero sign is the only thing the concordance statistic observes. -/
noncomputable def slopeOfSigns {n : ℕ} (sgn : Fin n → Bool) : Fin n → ℝ :=
  fun i => if sgn i then 1 else -1

/-- The number of trajectories whose sign is positive. -/
def positiveCount {n : ℕ} (sgn : Fin n → Bool) : ℕ :=
  (Finset.univ.filter fun i => sgn i = true).card

/-- The number of trajectories whose sign is negative. -/
def negativeCount {n : ℕ} (sgn : Fin n → Bool) : ℕ :=
  (Finset.univ.filter fun i => sgn i = false).card

theorem slopeOfSigns_ne_zero {n : ℕ} (sgn : Fin n → Bool) (i : Fin n) :
    slopeOfSigns sgn i ≠ 0 := by
  unfold slopeOfSigns
  cases sgn i <;> simp

theorem directionSign_slopeOfSigns {n : ℕ} (sgn : Fin n → Bool) (i : Fin n) :
    directionSign (slopeOfSigns sgn i) = (if sgn i then 1 else -1) := by
  cases h : sgn i <;> simp [slopeOfSigns, directionSign, h]

theorem filter_directionSign_eq_one {n : ℕ} (sgn : Fin n → Bool) :
    (Finset.univ.filter fun i => directionSign (slopeOfSigns sgn i) = 1) =
      (Finset.univ.filter fun i => sgn i = true) := by
  apply Finset.filter_congr
  intro i _
  rw [directionSign_slopeOfSigns]
  cases sgn i <;> simp

theorem filter_directionSign_eq_neg_one {n : ℕ} (sgn : Fin n → Bool) :
    (Finset.univ.filter fun i => directionSign (slopeOfSigns sgn i) = -1) =
      (Finset.univ.filter fun i => sgn i = false) := by
  apply Finset.filter_congr
  intro i _
  rw [directionSign_slopeOfSigns]
  cases sgn i <;> simp

theorem positiveCount_add_negativeCount {n : ℕ} (sgn : Fin n → Bool) :
    positiveCount sgn + negativeCount sgn = n := by
  have hf : (Finset.univ.filter fun i => sgn i = false) =
      (Finset.univ.filter fun a => ¬ sgn a = true) := by
    apply Finset.filter_congr
    intro i _
    cases sgn i <;> simp
  have h := Finset.card_filter_add_card_filter_not
    (s := (Finset.univ : Finset (Fin n))) (p := fun i => sgn i = true)
  simp only [Finset.card_univ, Fintype.card_fin] at h
  unfold positiveCount negativeCount
  rw [hf, h]

/-- Under a sign vector, the directional concordance depends only on the number of
    positive signs: it is the fair folded-majority `fairNullConcordance`. -/
theorem directionalConcordance_slopeOfSigns {n : ℕ} (hn : 0 < n) (sgn : Fin n → Bool) :
    directionalConcordance (Finset.univ : Finset (Fin n)) (slopeOfSigns sgn) =
      some (fairNullConcordance n (positiveCount sgn)) := by
  rw [directionalConcordance_eq_counts]
  simp only [concordanceCounts]
  rw [filter_directionSign_eq_one, filter_directionSign_eq_neg_one]
  change (if positiveCount sgn + negativeCount sgn = 0 then none
      else some (max (positiveCount sgn : ℝ) (negativeCount sgn : ℝ) /
        ((positiveCount sgn + negativeCount sgn : ℕ) : ℝ))) =
      some (fairNullConcordance n (positiveCount sgn))
  have hsum : positiveCount sgn + negativeCount sgn = n := positiveCount_add_negativeCount sgn
  have hneg : negativeCount sgn = n - positiveCount sgn := by omega
  rw [hsum, if_neg (Nat.ne_of_gt hn), hneg]
  have hp : positiveCount sgn ≤ n := by omega
  simp only [fairNullConcordance, Nat.cast_max, Nat.cast_sub hp]

/-- The concordance statistic is a *symmetric* function of the signs: two sign vectors
    with the same positive count have identical concordance. This is the exchangeability
    consequence used to reduce the null law to a distribution on the count. -/
theorem directionalConcordance_eq_of_sameCount {n : ℕ} (hn : 0 < n)
    {sgn₁ sgn₂ : Fin n → Bool} (h : positiveCount sgn₁ = positiveCount sgn₂) :
    directionalConcordance (Finset.univ : Finset (Fin n)) (slopeOfSigns sgn₁) =
      directionalConcordance (Finset.univ : Finset (Fin n)) (slopeOfSigns sgn₂) := by
  rw [directionalConcordance_slopeOfSigns hn sgn₁,
    directionalConcordance_slopeOfSigns hn sgn₂, h]

/-! ## Law of total probability for an arbitrary sign weight -/

/-- Expected concordance under an arbitrary finite sign weight. The value is the mean of
    `fairNullConcordance n (positiveCount ·)`, the concordance of the sign vector. -/
noncomputable def expectedConcordance (n : ℕ) (weight : (Fin n → Bool) → ℝ) : ℝ :=
  ∑ sgn : Fin n → Bool, weight sgn * fairNullConcordance n (positiveCount sgn)

/-- The sign-vector space is canonically a finite set of *positive-sign subsets*. -/
def signEquivSubset (n : ℕ) : (Fin n → Bool) ≃ Finset (Fin n) where
  toFun sgn := Finset.univ.filter fun i => sgn i = true
  invFun S := fun i => decide (i ∈ S)
  left_inv sgn := by
    funext i
    cases h : sgn i <;> simp [h]
  right_inv S := by
    ext i
    simp

theorem positiveCount_eq_card (n : ℕ) (sgn : Fin n → Bool) :
    positiveCount sgn = (signEquivSubset n sgn).card := rfl

/-- **Law of total probability (regrouping).** Summing any weight over sign vectors
    groups into a sum over the subsets of positive signs, because a sign vector is a
    subset in disguise. No exchangeability assumption is needed for this identity. -/
theorem expectedConcordance_eq_subsetSum (n : ℕ) (weight : (Fin n → Bool) → ℝ)
    (w : Finset (Fin n) → ℝ) (hw : ∀ sgn, weight sgn = w (signEquivSubset n sgn)) :
    expectedConcordance n weight =
      ∑ S : Finset (Fin n), w S * fairNullConcordance n S.card := by
  unfold expectedConcordance
  have h : ∑ sgn : Fin n → Bool, weight sgn * fairNullConcordance n (positiveCount sgn) =
      ∑ sgn : Fin n → Bool,
        (fun S => w S * fairNullConcordance n S.card) (signEquivSubset n sgn) := by
    apply Finset.sum_congr rfl
    intro sgn _
    rw [hw sgn, positiveCount_eq_card]
  rw [h]
  exact Fintype.sum_equiv (signEquivSubset n) _ _ (fun _ => rfl)

/-- The uniform (exchangeable, iid fair-sign) null expectation. -/
noncomputable def uniformNullExpectation (n : ℕ) : ℝ :=
  expectedConcordance n (fun _ => (1 / 2 : ℝ) ^ n)

/-- **Uniform null equals the fair folded-binomial mean.** The uniform sign measure is the
    exchangeable null; conditioning on the positive count and counting subsets of each size
    gives `fairNullMean`. This is the law-of-total-probability content of the null. -/
theorem uniformNullExpectation_eq_fairNullMean (n : ℕ) :
    uniformNullExpectation n = fairNullMean n := by
  unfold uniformNullExpectation
  rw [expectedConcordance_eq_subsetSum n _ (fun _ => (1 / 2 : ℝ) ^ n) (fun _ => rfl)]
  rw [show (∑ S : Finset (Fin n), (1 / 2 : ℝ) ^ n * fairNullConcordance n S.card)
      = ∑ S ∈ (Finset.univ : Finset (Fin n)).powerset,
          (1 / 2 : ℝ) ^ n * fairNullConcordance n S.card by
        rw [Finset.powerset_univ]]
  have h := Finset.sum_powerset_apply_card
    (f := fun k => (1 / 2 : ℝ) ^ n * fairNullConcordance n k)
    (x := (Finset.univ : Finset (Fin n)))
  rw [h]
  simp only [Finset.card_univ, Fintype.card_fin]
  rw [fairNullMean, fairBinomialAverage]
  apply Finset.sum_congr rfl
  intro k _
  rw [nsmul_eq_mul, fairBinomialWeight, div_pow, one_pow]
  ring

end TrajectorySelection
