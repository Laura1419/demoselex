import Mathlib.Tactic

namespace TrajectorySelection
namespace Collinearity

/-!
# Time–GRM alignment in a linearized mixed model

The time and eigenvector coordinates below are **already residualized against the
intercept and any other fixed covariates**. A unit vector `u` represents one
GRM eigenmode, with its sign immaterial. `residual t u` is the hard projection
of time off that mode; this is a limiting/high-variance-mode idealization, not
the exact estimator of a finite-variance binomial GLMM.
-/

/-- Euclidean dot product of finite observation vectors. -/
def dot {n : ℕ} (a b : Fin n → ℝ) : ℝ := ∑ i, a i * b i

/-- Squared Euclidean length (the unweighted time sum of squares). -/
def energy {n : ℕ} (a : Fin n → ℝ) : ℝ := dot a a

/-- Squared length cannot be negative (also for the projected time vector). -/
theorem energy_nonneg {n : ℕ} (t : Fin n → ℝ) : 0 ≤ energy t := by
  unfold energy dot
  exact Finset.sum_nonneg fun i _ => mul_self_nonneg (t i)

/-- Time remaining after removing the component along a *unit* GRM eigenvector. -/
def residual {n : ℕ} (t u : Fin n → ℝ) : Fin n → ℝ :=
  fun i => t i - dot t u * u i

/-- Squared cosine with a unit eigenvector; defined for nonzero time energy. -/
noncomputable def cosineSquared {n : ℕ} (t u : Fin n → ℝ) : ℝ :=
  dot t u ^ 2 / energy t

/-- The information-variance ratio under hard removal of the mode.
    Only interpreted when the residual energy is strictly positive. -/
noncomputable def projectedVIF {n : ℕ} (t u : Fin n → ℝ) : ℝ :=
  energy t / energy (residual t u)

/-- PQL/GLS information for one random-effect eigenmode in the simplified
    working covariance `σ² I + λ uuᵀ`, where `λ ≥ 0` is the GRM-mode
    variance and `σ² > 0` is the isotropic working-noise variance.
    Unlike a hard projection, finite `λ` retains the parallel component. -/
noncomputable def workingInformation {n : ℕ} (t u : Fin n → ℝ)
    (sigma2 modeVariance : ℝ) : ℝ :=
  energy (residual t u) / sigma2 + dot t u ^ 2 / (sigma2 + modeVariance)

/-- Pythagoras for a unit eigenvector: the lost time sum of squares is
    exactly its squared loading on that eigenvector. -/
theorem residual_energy {n : ℕ} (t u : Fin n → ℝ) (hu : energy u = 1) :
    energy (residual t u) = energy t - dot t u ^ 2 := by
  have expand : energy (residual t u) =
      energy t - 2 * dot t u * dot t u + dot t u ^ 2 * energy u := by
    simp only [energy, residual, dot]
    calc
      (∑ i : Fin n, (t i - (∑ j : Fin n, t j * u j) * u i) *
          (t i - (∑ j : Fin n, t j * u j) * u i)) =
          ∑ i : Fin n, (t i * t i - 2 * (∑ j : Fin n, t j * u j) *
            (t i * u i) + (∑ j : Fin n, t j * u j) ^ 2 * (u i * u i)) := by
              apply Finset.sum_congr rfl
              intro i _
              ring
      _ = _ := by
        simp only [Finset.sum_sub_distrib, Finset.sum_add_distrib,
          ← Finset.mul_sum]
  rw [expand, hu]
  ring

/-- The within-mode residual fraction is `1 − cos² θ`. -/
theorem residual_fraction {n : ℕ} (t u : Fin n → ℝ)
    (hu : energy u = 1) (ht : energy t ≠ 0) :
    energy (residual t u) / energy t = 1 - cosineSquared t u := by
  rw [residual_energy t u hu]
  unfold cosineSquared
  rw [sub_div]
  simp [ht]

/-- Equivalently, the projected time sum of squares collapses to
    `(1 − cos² θ)` times the original sum of squares. -/
theorem residual_energy_cosine {n : ℕ} (t u : Fin n → ℝ)
    (hu : energy u = 1) (ht : energy t ≠ 0) :
    energy (residual t u) = energy t * (1 - cosineSquared t u) := by
  have h := residual_fraction t u hu ht
  apply (div_eq_iff ht).mp at h
  linear_combination h

/-- The squared cosine really lies in the unit interval when time has positive energy. -/
theorem cosineSquared_nonneg {n : ℕ} (t u : Fin n → ℝ)
    (ht : 0 < energy t) : 0 ≤ cosineSquared t u := by
  unfold cosineSquared
  exact div_nonneg (sq_nonneg _) (le_of_lt ht)

theorem cosineSquared_le_one {n : ℕ} (t u : Fin n → ℝ)
    (hu : energy u = 1) (ht : 0 < energy t) : cosineSquared t u ≤ 1 := by
  have h := residual_fraction t u hu (ne_of_gt ht)
  have hnon : 0 ≤ energy (residual t u) / energy t :=
    div_nonneg (energy_nonneg _) (le_of_lt ht)
  rw [h] at hnon
  linarith

/-- A near-parallel time vector leaves at most `ε` of its original time
    sum of squares after hard projection. -/
theorem near_alignment {n : ℕ} (t u : Fin n → ℝ) (epsilon : ℝ)
    (hu : energy u = 1) (ht : 0 < energy t)
    (hnear : 1 - cosineSquared t u ≤ epsilon) :
    energy (residual t u) ≤ epsilon * energy t := by
  rw [residual_energy_cosine t u hu (ne_of_gt ht)]
  nlinarith

/-- When the projected mode does not exhaust time, the hard-projection
    variance-inflation factor is `1 / (1 − cos² θ)`. -/
theorem projectedVIF_eq {n : ℕ} (t u : Fin n → ℝ)
    (hu : energy u = 1) (ht : energy t ≠ 0)
    (hr : energy (residual t u) ≠ 0) :
    projectedVIF t u = 1 / (1 - cosineSquared t u) := by
  have h := residual_energy_cosine t u hu ht
  have hc : 1 - cosineSquared t u ≠ 0 := by
    intro hz
    rw [hz, mul_zero] at h
    exact hr h
  unfold projectedVIF
  rw [h]
  field_simp

/-- Perfect collinearity has zero *projected* time information. This does
    not assert that a finite-variance mixed model is unidentifiable. -/
theorem perfect_alignment {n : ℕ} (t u : Fin n → ℝ)
    (hu : energy u = 1) (ht : energy t ≠ 0)
    (hcos : cosineSquared t u = 1) :
    energy (residual t u) = 0 := by
  rw [residual_energy_cosine t u hu ht, hcos]
  ring

/-- In the finite-variance working model, the mode's contribution is
    retained rather than deleted (and vanishes only in a high-variance limit). -/
theorem workingInformation_excess {n : ℕ} (t u : Fin n → ℝ)
    (sigma2 modeVariance : ℝ) :
    workingInformation t u sigma2 modeVariance - energy (residual t u) / sigma2 =
      dot t u ^ 2 / (sigma2 + modeVariance) := by
  unfold workingInformation
  ring

/-- The **finite** PQL/GLS information fraction relative to isotropic noise:
    the parallel fraction is attenuated by `σ²/(σ²+λ)`, not removed. -/
theorem working_information_fraction {n : ℕ} (t u : Fin n → ℝ)
    (sigma2 modeVariance : ℝ) (hu : energy u = 1)
    (ht : energy t ≠ 0) (hσ : sigma2 ≠ 0)
    (hden : sigma2 + modeVariance ≠ 0) :
    workingInformation t u sigma2 modeVariance / (energy t / sigma2) =
      (1 - cosineSquared t u) +
        cosineSquared t u * sigma2 / (sigma2 + modeVariance) := by
  have hr := residual_energy_cosine t u hu ht
  have hd : energy t * cosineSquared t u = dot t u ^ 2 := by
    unfold cosineSquared
    field_simp
  unfold workingInformation
  rw [hr, ← hd]
  field_simp

/-- If the remaining time fraction is at most `ε` and the retained fraction
    of the GRM mode is at most `δ`, finite-model information is at most
    `ε + δ` of its isotropic baseline. No finite mode is silently deleted. -/
theorem near_alignment_working_bound {n : ℕ} (t u : Fin n → ℝ)
    (sigma2 modeVariance epsilon delta : ℝ)
    (hu : energy u = 1) (ht : 0 < energy t)
    (hσ : 0 < sigma2) (hMode : 0 ≤ modeVariance)
    (hnear : 1 - cosineSquared t u ≤ epsilon)
    (hlarge : sigma2 / (sigma2 + modeVariance) ≤ delta) :
    workingInformation t u sigma2 modeVariance / (energy t / sigma2) ≤
      epsilon + delta := by
  have hden : 0 < sigma2 + modeVariance := add_pos_of_pos_of_nonneg hσ hMode
  have hcos := cosineSquared_le_one t u hu ht
  have hratio : 0 ≤ sigma2 / (sigma2 + modeVariance) := le_of_lt (div_pos hσ hden)
  have hproduct : cosineSquared t u * (sigma2 / (sigma2 + modeVariance)) ≤
      sigma2 / (sigma2 + modeVariance) := by
    nlinarith [mul_nonneg (sub_nonneg.mpr hcos) hratio]
  rw [working_information_fraction t u sigma2 modeVariance hu
    (ne_of_gt ht) (ne_of_gt hσ) (ne_of_gt hden)]
  rw [mul_div_assoc]
  linarith

/-- Exact alignment still carries positive information at finite mode variance;
    the assertion of zero identifiability applies only to hard projection or
    the `λ → ∞` limit. -/
theorem perfect_alignment_working {n : ℕ} (t u : Fin n → ℝ)
    (sigma2 modeVariance : ℝ) (hu : energy u = 1)
    (ht : energy t ≠ 0) (hcos : cosineSquared t u = 1) :
    workingInformation t u sigma2 modeVariance =
      energy t / (sigma2 + modeVariance) := by
  have hr := perfect_alignment t u hu ht hcos
  have hd : energy t * cosineSquared t u = dot t u ^ 2 := by
    unfold cosineSquared
    field_simp
  unfold workingInformation
  rw [hr, zero_div, zero_add, ← hd, hcos, mul_one]

/-- If the covariance is positive definite, the finite-mode information is
    at least the hard-projection information. -/
theorem workingInformation_ge_projected {n : ℕ} (t u : Fin n → ℝ)
    (sigma2 modeVariance : ℝ) (hσ : 0 < sigma2) (hMode : 0 ≤ modeVariance) :
    energy (residual t u) / sigma2 ≤ workingInformation t u sigma2 modeVariance := by
  have hden : 0 < sigma2 + modeVariance := add_pos_of_pos_of_nonneg hσ hMode
  unfold workingInformation
  exact le_add_of_nonneg_right (div_nonneg (sq_nonneg _) (le_of_lt hden))

end Collinearity
end TrajectorySelection
