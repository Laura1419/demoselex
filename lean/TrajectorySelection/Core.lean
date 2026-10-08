/-
Core definitions and basic types for trajectory replication

AXIOMS policy: Every `sorry`, `axiom`, and `native_decide` must be documented in lean/AXIOMS.md
-/

import Mathlib.Data.Real.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Data.Nat.Cast.Basic

namespace TrajectorySelection

/- Basic types for population genetics -/
abbrev PopulationId := Nat
abbrev TimePoint := Nat
abbrev AlleleFrequency := ℝ
abbrev LocusId := Nat

/- A trajectory is a sequence of allele frequencies over time for a population -/
structure Trajectory where
  population : PopulationId
  locus : LocusId
  frequencies : List (TimePoint × AlleleFrequency)

/- Basic trajectory operations -/
def Trajectory.length (t : Trajectory) : Nat := t.frequencies.length

def Trajectory.get_frequency (t : Trajectory) (time : TimePoint) : Option AlleleFrequency :=
  t.frequencies.find? (λ (tp, _) => tp = time) |>.map Prod.snd

/- Direction of frequency change -/
inductive Direction where
  | increasing : Direction
  | decreasing : Direction
  | stable : Direction

-- Derive BEq for Direction
deriving instance BEq for Direction

-- Derive Inhabited for Trajectory
deriving instance Inhabited for Trajectory

end TrajectorySelection