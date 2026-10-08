import TrajectorySelection.Defs
import Mathlib.Tactic

namespace TrajectorySelection

/-! # A two-trajectory sign experiment

The constant source frequencies below do **not** by themselves imply a random-sign
null. The additional, explicitly calibrated observation model assumes that the two
estimated slope signs are independent fair Bernoulli variables under replacement.
Under replication, positive true drift is accompanied by an *assumed* sign-detection
probability `q > 1/2`, independently for both trajectories. The implication from a
physical drift/noise law to this probability is outside this finite sign model.
-/

/-- Scenario A: two temporally constant sources, with arbitrary time-dependent
mixture. `observedSigns` are signs of noisy *estimated* within-source slopes, not
signs of the true (zero) within-source changes. -/
structure ReplacementScenario (Time : Type*) where
  sources : ReplacementOnly Time
  observedSigns : Bool × Bool

/-- The true frequency along each source is fixed across all times. -/
def ReplacementScenario.withinSource {Time : Type*}
    (scenario : ReplacementScenario Time) (trajectory : Fin 2) (_time : Time) : ℝ :=
  if trajectory = 0 then scenario.sources.sourceA else scenario.sources.sourceB

theorem ReplacementScenario.withinSource_constant {Time : Type*}
    (scenario : ReplacementScenario Time) (trajectory : Fin 2) (t u : Time) :
    scenario.withinSource trajectory t = scenario.withinSource trajectory u := rfl

/-- Each source produces a bounded constant path; the pooled mixture need not be constant. -/
def ReplacementScenario.path {Time : Type*}
    (scenario : ReplacementScenario Time) (trajectory : Fin 2) :
    AlleleFrequencyPath Time where
  frequency := scenario.withinSource trajectory
  bounded := by
    intro time
    unfold ReplacementScenario.withinSource
    split_ifs
    · exact scenario.sources.sourceA_bounded
    · exact scenario.sources.sourceB_bounded

/-- Scenario B: two paths move by the same strictly positive amount between
`false` (early) and `true` (late); the bounds make these valid frequencies. -/
structure ReplicatedChange where
  initial : Fin 2 → ℝ
  drift : ℝ
  drift_pos : 0 < drift
  initial_nonneg : ∀ i, 0 ≤ initial i
  final_le_one : ∀ i, initial i + drift ≤ 1
  /-- Assumed probability that a noisy slope estimate reports the correct sign. -/
  detectionProbability : ℝ
  detection_biased : 1 / 2 < detectionProbability
  detection_le_one : detectionProbability ≤ 1

/-- Exact true path under replicated directional change. -/
def ReplicatedChange.frequency (scenario : ReplicatedChange)
    (trajectory : Fin 2) (late : Bool) : ℝ :=
  scenario.initial trajectory + if late then scenario.drift else 0

theorem ReplicatedChange.change (scenario : ReplicatedChange) (i : Fin 2) :
    scenario.frequency i true - scenario.frequency i false = scenario.drift := by
  simp [ReplicatedChange.frequency]

theorem ReplicatedChange.frequency_bounded (scenario : ReplicatedChange)
    (i : Fin 2) (late : Bool) :
    0 ≤ scenario.frequency i late ∧ scenario.frequency i late ≤ 1 := by
  cases late <;> simp only [ReplicatedChange.frequency, Bool.false_eq_true, ite_false,
    ite_true, add_zero]
  · constructor
    · exact scenario.initial_nonneg i
    · linarith [scenario.final_le_one i, scenario.drift_pos]
  · constructor
    · linarith [scenario.initial_nonneg i, scenario.drift_pos]
    · exact scenario.final_le_one i

/-- The observed bounded frequency path on one changing trajectory. -/
def ReplicatedChange.path (scenario : ReplicatedChange) (i : Fin 2) :
    AlleleFrequencyPath Bool where
  frequency := scenario.frequency i
  bounded := scenario.frequency_bounded i

/-- A sign pair interpreted as two nonzero estimated slopes. -/
def signSlope (signs : Bool × Bool) (i : Fin 2) : ℝ :=
  if i = 0 then (if signs.1 then 1 else -1) else (if signs.2 then 1 else -1)

/-- The existing concordance statistic, evaluated on both trajectory labels. -/
noncomputable def pairConcordance (signs : Bool × Bool) : Option ℝ :=
  directionalConcordance (Finset.univ : Finset (Fin 2)) (signSlope signs)

/-- The statistic observed in a particular replacement scenario. -/
noncomputable def ReplacementScenario.concordance {Time : Type*}
    (scenario : ReplacementScenario Time) : Option ℝ :=
  pairConcordance scenario.observedSigns

/-- At fixed observed signs the statistic cannot depend on source frequencies
or mixing weights. This is a claim about the sign-only statistic, not about how
real-world frequency estimates acquire their signs. -/
theorem replacement_concordance_determined_by_signs {Time : Type*}
    (a b : ReplacementOnly Time) (signs : Bool × Bool) :
    (ReplacementScenario.mk a signs).concordance =
      (ReplacementScenario.mk b signs).concordance := rfl

/-- The mass of two independent sign observations, each positive with probability `q`.
This is an algebraic finite probability model, not a claim about arbitrary errors. -/
def signPairMass (q : ℝ) (signs : Bool × Bool) : ℝ :=
  (if signs.1 then q else 1 - q) * (if signs.2 then q else 1 - q)

/-- The four outcome masses sum to one. Nonnegativity additionally requires `0 ≤ q ≤ 1`. -/
theorem signPairMass_total (q : ℝ) :
    (∑ signs : Bool × Bool, signPairMass q signs) = 1 := by
  simp only [Fintype.sum_prod_type, Fintype.sum_bool]
  simp [signPairMass]
  ring

/-- With a valid Bernoulli parameter these are nonnegative probability masses. -/
theorem signPairMass_nonneg (q : ℝ) (hq₀ : 0 ≤ q) (hq₁ : q ≤ 1)
    (signs : Bool × Bool) : 0 ≤ signPairMass q signs := by
  rcases signs with ⟨a, b⟩
  cases a <;> cases b <;> simp [signPairMass] <;> positivity

/-- Probability (finite sum of masses) of perfect two-source concordance. -/
noncomputable def highConcordanceMass (q : ℝ) : ℝ :=
  ∑ signs : Bool × Bool,
    if pairConcordance signs = some 1 then signPairMass q signs else 0

/-- Scenario A assigns iid fair signs, so its entire pushforward distribution of
concordance agrees pointwise with the exchangeable-sign null. Constancy of the
underlying paths alone is not used to manufacture this distribution: fairness and
independence are explicit in the mass `signPairMass (1/2)`. -/
theorem replacement_null_distribution {Time : Type*}
    (sources : ReplacementOnly Time) (event : Option ℝ → Prop) [DecidablePred event] :
    (∑ signs : Bool × Bool,
      if event ((ReplacementScenario.mk sources signs).concordance) then
        signPairMass (1 / 2) signs else 0) =
    (∑ signs : Bool × Bool,
      if event (pairConcordance signs) then (1 / 4 : ℝ) else 0) := by
  apply Finset.sum_congr rfl
  intro signs _
  change (if event (pairConcordance signs) then signPairMass (1 / 2) signs else 0) = _
  have h : signPairMass (1 / 2) signs = (1 / 4 : ℝ) := by
    cases signs with
    | mk x y => cases x <;> cases y <;> norm_num [signPairMass]
  by_cases he : event (pairConcordance signs)
  · simpa only [he, ↓reduceIte, one_div] using h
  · simp [he]

/-- Exact event mass: either both observed signs are positive, or both negative. -/
theorem highConcordanceMass_eq (q : ℝ) :
    highConcordanceMass q = q * q + (1 - q) * (1 - q) := by
  -- The four possible sign pairs are finite; mixed signs have concordance 1/2.
  have hu : (Finset.univ : Finset (Fin 2)) = {0, 1} := by decide
  simp only [highConcordanceMass, Fintype.sum_prod_type, Fintype.sum_bool]
  norm_num [pairConcordance, directionalConcordance, hu, signSlope,
    directionSign, signPairMass, Finset.filter_insert, Finset.filter_singleton]

/-- The fair-sign replacement null has perfect concordance with probability 1/2. -/
theorem replacement_high_mass : highConcordanceMass (1 / 2) = (1 / 2 : ℝ) := by
  rw [highConcordanceMass_eq]
  norm_num

/-- Scenario B: under independent sign detections with a common `q > 1/2`,
the probability of perfect concordance strictly exceeds the replacement null.
The assumed link between positive drift and `q > 1/2` is indispensable. -/
theorem replication_strictly_exceeds_null (scenario : ReplicatedChange) :
    highConcordanceMass scenario.detectionProbability >
      highConcordanceMass (1 / 2) := by
  rw [highConcordanceMass_eq, replacement_high_mass]
  nlinarith [scenario.detection_biased]

end TrajectorySelection
