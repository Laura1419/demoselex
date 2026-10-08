import Lake
open Lake DSL

-- Repository-root entry point; the canonical Lean project is in lean/.
package trajectory_replication

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.34.1"

@[default_target]
lean_lib TrajectorySelection where
  srcDir := "lean"
