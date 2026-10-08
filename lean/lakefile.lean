import Lake
open Lake DSL

-- Canonical project for the formalization in this directory.
package trajectory_replication

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.34.1"

@[default_target]
lean_lib TrajectorySelection
