import TrajectorySelection.Defs
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Tactic

namespace TrajectorySelection

/-- Concordance from the positive and negative *informative* trajectory counts.
    Undefined if there are no informative trajectories. -/
noncomputable def concordanceCounts (positive negative : ℕ) : Option ℝ :=
  if positive + negative = 0 then none
  else some ((max positive negative : ℝ) / (positive + negative : ℕ))

theorem directionalConcordance_eq_counts {T : Type*} (s : Finset T) (slope : T → ℝ) :
    directionalConcordance s slope =
      concordanceCounts (s.filter fun t => directionSign (slope t) = 1).card
        (s.filter fun t => directionSign (slope t) = -1).card := by
  simp only [directionalConcordance, concordanceCounts, Nat.cast_max]

/-- Majority concordance lies between one half and one, including a tie and unanimity. -/
theorem concordanceCounts_range (p q : ℕ) (h : 0 < p + q) :
    ∃ c, concordanceCounts p q = some c ∧ (1 / 2 : ℝ) ≤ c ∧ c ≤ 1 := by
  refine ⟨(max p q : ℝ) / (p + q : ℕ), ?_, ?_, ?_⟩
  · simp only [concordanceCounts, Nat.ne_of_gt h, ↓reduceIte]
  · have hn : (0 : ℝ) < (p + q : ℕ) := by exact_mod_cast h
    have hm : p + q ≤ 2 * max p q := by omega
    have hm' : (p : ℝ) + q ≤ 2 * (max p q : ℝ) := by exact_mod_cast hm
    have hs : (p : ℝ) + q = (↑(p + q) : ℝ) := by norm_cast
    rw [hs] at hm'
    apply (le_div_iff₀ hn).2
    nlinarith
  · have hn : (0 : ℝ) < (p + q : ℕ) := by exact_mod_cast h
    apply (div_le_iff₀ hn).2
    simpa only [one_mul] using
      (show (max p q : ℝ) ≤ (p + q : ℕ) by exact_mod_cast (show max p q ≤ p + q by omega))

/-- At fixed informative count, making the majority larger cannot lower concordance. -/
theorem concordance_majority_mono (n a b : ℕ) (hn : 0 < n) (hab : a ≤ b) :
    (a : ℝ) / n ≤ (b : ℝ) / n := by
  exact div_le_div_of_nonneg_right (by exact_mod_cast hab) (by exact_mod_cast hn.le)

/-- Correlation mass of distinct trajectory *families*. The diagonal is one; off-diagonal
    relatedness discounts the number of independent replicates. -/
def correlationMass {Family : Type*} (groups : Finset Family) (ρ : Family → Family → ℝ) : ℝ :=
  ∑ i ∈ groups, ∑ j ∈ groups, ρ i j

/-- Effective count of distinct trajectory families, `(number of families)^2 / correlation mass`.
    The zero-family case is undefined, rather than dividing by zero. `ρ` is assumed to have
    unit diagonal and values in `[0,1]` in the theorems below. -/
noncomputable def effectiveFamilyCount {Family : Type*}
    (groups : Finset Family) (ρ : Family → Family → ℝ) : Option ℝ :=
  if groups.card = 0 then none
  else some ((groups.card : ℝ) ^ 2 / correlationMass groups ρ)

/-- Assign every sampled trajectory a genealogical/relatedness family. Near-identical copies
    are assigned to the *same* family before the covariance discount is computed. This
    deduplication must be supplied by the inference model; numerical similarity alone does
    not imply a transitive clustering. -/
noncomputable def effectiveTrajectoryCount {T Family : Type*} [DecidableEq Family]
    (selected : Finset T) (family : T → Family) (ρ : Family → Family → ℝ) : Option ℝ :=
  effectiveFamilyCount (selected.image family) ρ

/-- An added trajectory assigned to an existing family cannot increase effective count:
    it leaves it exactly unchanged, even if its individual label is new. -/
theorem effectiveTrajectoryCount_duplicate {T Family : Type*} [DecidableEq Family]
    [DecidableEq T] (selected : Finset T) (family : T → Family)
    (ρ : Family → Family → ℝ) (copy : T)
    (hcopy : family copy ∈ selected.image family) :
    effectiveTrajectoryCount (insert copy selected) family ρ =
      effectiveTrajectoryCount selected family ρ := by
  simp [effectiveTrajectoryCount, Finset.image_insert, hcopy]

/-- With unit diagonal and correlations in `[0,1]`, total mass is between
    `m` and `m²`, where `m` is the number of distinct families. -/
theorem correlationMass_bounds {Family : Type*} (groups : Finset Family)
    (ρ : Family → Family → ℝ) (hdiag : ∀ i ∈ groups, ρ i i = 1)
    (hrange : ∀ i ∈ groups, ∀ j ∈ groups, 0 ≤ ρ i j ∧ ρ i j ≤ 1) :
    (groups.card : ℝ) ≤ correlationMass groups ρ ∧
      correlationMass groups ρ ≤ (groups.card : ℝ) ^ 2 := by
  classical
  have hlo (i : Family) (hi : i ∈ groups) :
      1 ≤ ∑ j ∈ groups, ρ i j := by
    calc
      (1 : ℝ) = ∑ j ∈ groups, if j = i then (1 : ℝ) else 0 := by simp [hi]
      _ ≤ ∑ j ∈ groups, ρ i j := by
        apply Finset.sum_le_sum
        intro j hj
        by_cases hij : j = i
        · simp [hij, hdiag i hi]
        · simpa [hij] using (hrange i hi j hj).1
  have hhi (i : Family) (hi : i ∈ groups) :
      (groups.card : ℝ) ≥ ∑ j ∈ groups, ρ i j := by
    calc
      (groups.card : ℝ) = ∑ _j ∈ groups, (1 : ℝ) := by simp
      _ ≥ ∑ j ∈ groups, ρ i j := by
        apply Finset.sum_le_sum
        intro j hj
        exact (hrange i hi j hj).2
  constructor
  · calc
      (groups.card : ℝ) = ∑ _i ∈ groups, (1 : ℝ) := by simp
      _ ≤ correlationMass groups ρ := by
        unfold correlationMass
        apply Finset.sum_le_sum
        exact hlo
  · calc
      correlationMass groups ρ ≤ ∑ _i ∈ groups, (groups.card : ℝ) := by
        unfold correlationMass
        apply Finset.sum_le_sum
        exact hhi
      _ = (groups.card : ℝ) ^ 2 := by simp [pow_two]

/-- For a nonempty family set the effective count lies in `[1,m]`. -/
theorem effectiveFamilyCount_range {Family : Type*} (groups : Finset Family)
    (ρ : Family → Family → ℝ) (hn : 0 < groups.card)
    (hdiag : ∀ i ∈ groups, ρ i i = 1)
    (hrange : ∀ i ∈ groups, ∀ j ∈ groups, 0 ≤ ρ i j ∧ ρ i j ≤ 1) :
    ∃ ne, effectiveFamilyCount groups ρ = some ne ∧
      1 ≤ ne ∧ ne ≤ (groups.card : ℝ) := by
  have hm := correlationMass_bounds groups ρ hdiag hrange
  have hcard : (0 : ℝ) < groups.card := by exact_mod_cast hn
  have hmass : 0 < correlationMass groups ρ := lt_of_lt_of_le hcard hm.1
  refine ⟨(groups.card : ℝ) ^ 2 / correlationMass groups ρ, ?_, ?_, ?_⟩
  · simp [effectiveFamilyCount, Nat.ne_of_gt hn]
  · apply (le_div_iff₀ hmass).2
    simpa only [one_mul] using hm.2
  · apply (div_le_iff₀ hmass).2
    have h : (0 : ℝ) ≤ (groups.card : ℝ) := hcard.le
    nlinarith [mul_nonneg h (sub_nonneg.mpr hm.1)]

/-- Increasing all pairwise relatedness on the same nonempty family set cannot increase
    the effective number. Positive masses are essential for this antitonicity. -/
theorem effectiveFamilyCount_antitone {Family : Type*} (groups : Finset Family)
    (ρ σ : Family → Family → ℝ) (hn : 0 < groups.card)
    (hρ : 0 < correlationMass groups ρ)
    (hσ : 0 < correlationMass groups σ)
    (hmore : ∀ i ∈ groups, ∀ j ∈ groups, ρ i j ≤ σ i j) :
    ∃ a b, effectiveFamilyCount groups ρ = some a ∧
      effectiveFamilyCount groups σ = some b ∧ b ≤ a := by
  have hm : correlationMass groups ρ ≤ correlationMass groups σ := by
    unfold correlationMass
    apply Finset.sum_le_sum
    intro i hi
    apply Finset.sum_le_sum
    intro j hj
    exact hmore i hi j hj
  refine ⟨(groups.card : ℝ) ^ 2 / correlationMass groups ρ,
    (groups.card : ℝ) ^ 2 / correlationMass groups σ, ?_, ?_, ?_⟩
  · simp [effectiveFamilyCount, Nat.ne_of_gt hn]
  · simp [effectiveFamilyCount, Nat.ne_of_gt hn]
  · apply (div_le_div_iff₀ hσ hρ).2
    have hc : 0 ≤ (groups.card : ℝ) ^ 2 := sq_nonneg _
    nlinarith [mul_nonneg hc (sub_nonneg.mpr hm)]

/-- Conditional on exactly `n` informative trajectories, the null with *jointly uniform*
    direction vectors (independent fair signs) has `K ∼ Binomial(n, 1/2)`. Exchangeability
    by itself does not imply this probability mass function. -/
noncomputable def fairBinomialWeight (n k : ℕ) : ℝ :=
  (n.choose k : ℝ) / 2 ^ n

noncomputable def fairBinomialAverage (n : ℕ) (f : ℕ → ℝ) : ℝ :=
  ∑ k ∈ Finset.range (n + 1), fairBinomialWeight n k * f k

theorem fairBinomialWeight_total (n : ℕ) :
    fairBinomialAverage n (fun _ => 1) = 1 := by
  simp only [fairBinomialAverage, mul_one, fairBinomialWeight, ← Finset.sum_div]
  rw [← Nat.cast_sum, Nat.sum_range_choose, Nat.cast_pow, Nat.cast_ofNat]
  exact div_self (pow_ne_zero _ (by norm_num))

/-- The positive count `k` has negative count `n-k`. Outside the support `k≤n`,
    this formula is not a null outcome and no claim is made about it. -/
noncomputable def fairNullConcordance (n k : ℕ) : ℝ :=
  (max k (n - k) : ℝ) / n

/-- The majority statistic is the folded binomial count, not `K/n` (whose mean
    would be one half by symmetry). -/
theorem fairNullConcordance_folded (n k : ℕ) (hn : 0 < n) (hk : k ≤ n) :
    fairNullConcordance n k = (1 / 2 : ℝ) + |(k : ℝ) - n / 2| / n := by
  have hsub : (↑(n - k) : ℝ) = (n : ℝ) - k := Nat.cast_sub hk
  have hnr : (0 : ℝ) < n := by exact_mod_cast hn
  by_cases h : 2 * k ≤ n
  · have hr : k ≤ n - k := by omega
    have habs : (k : ℝ) - n / 2 ≤ 0 := by
      have hh : (2 : ℝ) * k ≤ n := by exact_mod_cast h
      linarith
    have hr' : (k : ℝ) ≤ (n : ℝ) - k := by
      rw [← hsub]
      exact_mod_cast hr
    unfold fairNullConcordance
    rw [max_eq_right hr', abs_of_nonpos habs]
    field_simp
    ring
  · have hr : n - k ≤ k := by omega
    have habs : 0 ≤ (k : ℝ) - n / 2 := by
      have hh : (n : ℝ) < 2 * k := by exact_mod_cast (show n < 2 * k by omega)
      linarith
    have hr' : (n : ℝ) - k ≤ k := by
      rw [← hsub]
      exact_mod_cast hr
    unfold fairNullConcordance
    rw [max_eq_left hr', abs_of_nonneg habs]
    field_simp
    ring

noncomputable def fairNullMean (n : ℕ) : ℝ :=
  fairBinomialAverage n (fairNullConcordance n)

noncomputable def fairNullVariance (n : ℕ) : ℝ :=
  fairBinomialAverage n (fun k => (fairNullConcordance n k - fairNullMean n) ^ 2)

/-- Exact finite-sum null mean. The binomial coefficients include both tails and the tie. -/
theorem fairNullMean_formula (n : ℕ) :
    fairNullMean n =
      (∑ k ∈ Finset.range (n + 1), (n.choose k : ℝ) * fairNullConcordance n k) /
        2 ^ n := by
  simp only [fairNullMean, fairBinomialAverage, fairBinomialWeight, div_mul_eq_mul_div,
    ← Finset.sum_div]

/-- The null mean has an upward majority-selection bias of the expected folded
    deviation. In particular it is generally *not* one half under fair signs. -/
theorem fairNullMean_folded (n : ℕ) (hn : 0 < n) :
    fairNullMean n = (1 / 2 : ℝ) +
      fairBinomialAverage n (fun k => |(k : ℝ) - n / 2| / n) := by
  have hnorm := fairBinomialWeight_total n
  calc
    fairNullMean n = fairBinomialAverage n
      (fun k => (1 / 2 : ℝ) + |(k : ℝ) - n / 2| / n) := by
        unfold fairNullMean fairBinomialAverage
        apply Finset.sum_congr rfl
        intro k hk
        rw [fairNullConcordance_folded n k hn (Nat.le_of_lt_succ (Finset.mem_range.mp hk))]
    _ = (1 / 2 : ℝ) * fairBinomialAverage n (fun _ => 1) +
        fairBinomialAverage n (fun k => |(k : ℝ) - n / 2| / n) := by
        simp only [fairBinomialAverage, mul_add, Finset.sum_add_distrib,
          Finset.mul_sum, mul_one]
        congr 1
        apply Finset.sum_congr rfl
        intro k _
        ring
    _ = _ := by rw [hnorm]; ring

/-- The null variance is the second raw moment minus the square of the exact mean. -/
theorem fairNullVariance_formula (n : ℕ) :
    fairNullVariance n =
      fairBinomialAverage n (fun k => fairNullConcordance n k ^ 2) -
        fairNullMean n ^ 2 := by
  let μ := fairNullMean n
  have hnorm := fairBinomialWeight_total n
  have hmean : fairBinomialAverage n (fairNullConcordance n) = μ := rfl
  change (∑ k ∈ Finset.range (n + 1),
      fairBinomialWeight n k * (fairNullConcordance n k - μ) ^ 2) = _
  change (∑ k ∈ Finset.range (n + 1), fairBinomialWeight n k * 1) = 1 at hnorm
  change (∑ k ∈ Finset.range (n + 1),
      fairBinomialWeight n k * fairNullConcordance n k) = μ at hmean
  calc
    ∑ k ∈ Finset.range (n + 1), fairBinomialWeight n k * (fairNullConcordance n k - μ) ^ 2 =
      (∑ k ∈ Finset.range (n + 1), fairBinomialWeight n k * fairNullConcordance n k ^ 2) -
        2 * μ * (∑ k ∈ Finset.range (n + 1), fairBinomialWeight n k * fairNullConcordance n k) +
          μ ^ 2 * (∑ k ∈ Finset.range (n + 1), fairBinomialWeight n k * 1) := by
        rw [Finset.mul_sum, Finset.mul_sum]
        rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
        apply Finset.sum_congr rfl
        intro k _
        ring
    _ = _ := by
      rw [hmean, hnorm]
      change _ = (∑ k ∈ Finset.range (n + 1),
        fairBinomialWeight n k * fairNullConcordance n k ^ 2) - μ ^ 2
      ring

end TrajectorySelection
