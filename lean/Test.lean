/-
Basic test file to verify the Lean project builds correctly
-/

import TrajectorySelection

open TrajectorySelection

namespace Test

/- Create a simple test trajectory -/
def test_trajectory : Trajectory := {
  population := 1
  locus := 1001
  frequencies := [(0, 0.1), (10, 0.3), (20, 0.5)]
}

/- Test basic operations -/
noncomputable def test_direction : Option Direction := direction_of_change test_trajectory 0 20

noncomputable def test_monotonic : Bool :=
  let fp := FrequencyPath.fromTrajectory test_trajectory
  fp.isMonotonic

end Test