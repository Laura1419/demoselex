/-
Trajectory operations and analysis

AXIOMS policy: Every `sorry`, `axiom`, and `native_decide` must be documented in lean/AXIOMS.md
-/

import TrajectorySelection.Core

namespace TrajectorySelection

/- Calculate direction of frequency change between two time points -/
noncomputable def direction_of_change (t : Trajectory) (t1 t2 : TimePoint) : Option Direction := do
  let f1 ← t.get_frequency t1
  let f2 ← t.get_frequency t2
  if f1 < f2 then
    pure Direction.increasing
  else if f1 > f2 then
    pure Direction.decreasing
  else
    pure Direction.stable

/- Check if trajectory shows consistent direction over multiple time intervals -/
noncomputable def consistently_directed (t : Trajectory) (time_pairs : List (TimePoint × TimePoint)) : Bool :=
  let directions := time_pairs.filterMap (λ (t1, t2) => direction_of_change t t1 t2)
  let unique_directions := directions.eraseDups
  unique_directions.length ≤ 1

/- Calculate rate of frequency change -/
noncomputable def rate_of_change (t : Trajectory) (t1 t2 : TimePoint) : Option ℝ := do
  let f1 ← t.get_frequency t1
  let f2 ← t.get_frequency t2
  let time_diff : ℝ := ↑(t2 - t1)  -- Use explicit cast via `↑`
  if time_diff > 0 then
    pure ((f2 - f1) / time_diff)
  else
    none

end TrajectorySelection