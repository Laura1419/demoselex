/-
Lean 4 project with Mathlib for trajectory replication research
-/

import Lake
open Lake DSL

package "trajectory-replication" {
  -- add package configuration options here
}

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "master"

@[default_target]
lean_lib "TrajectorySelection" {
  -- Library files are in the lean/ subdirectory
  srcDir := "lean"
}