# Lean trust and gap ledger

## Accounting policy

For every occurrence of Lean `sorry`, declaration of `axiom`, or tactic `native_decide` in `lean/**/*.lean`, add a ledger entry **in the same commit**. Include the file and declaration, exact proposition/role, reason it is needed, assumptions or external computation trusted, impact on downstream claims, owner/status, and a concrete discharge or replacement plan. No unlisted use is permitted, including temporary workarounds. Review the ledger against a source search before accepting a proof; inspect `#print axioms <theorem>` for important results to identify inherited axioms as well as local gaps. Mathlib's own trusted foundations are an external dependency, not a license to add project axioms silently.

`sorry` means the statement is *not proved*, even though Lean can compile it. An `axiom` makes an unproved assumption part of the theory. `native_decide` delegates a decision to native evaluation rather than kernel reduction and must be called out as such. Neither a successful `lake build` nor a theorem that depends on these should be reported as a completed mathematical proof. Document changes and their scientific implications here, and remove entries only when the corresponding use is gone (retain a note in the journal).

## Inventory (initial bootstrap)

| File / declaration | Construct | Purpose / assumptions / impact | Status and discharge |
| --- | --- | --- | --- |
| None | `sorry`, `axiom`, `native_decide` | No project uses at bootstrap. | Re-audit on every Lean change. |

The following are **not** Lean axioms or `sorry`s, but are scientific placeholders in the inherited scaffold and must not be interpreted as calibrated statistics or validated formal definitions:

| Location | Current behavior | Necessary replacement |
| --- | --- | --- |
| `TrajectorySelection/FrequencyPath.lean`, `concordance_score` | Returns constant `some 0.5`, ignoring both inputs. | Specify a concordance statistic, handle missing/non-aligned time points, and prove its intended properties. |
| `TrajectorySelection/Concordance.lean`, `check_replication` | Compares a count of non-stable directions with a provisional `0.7` cutoff (not actual sign concordance). | Define replication over distinct histories and calibrate it under a null model. |
| `TrajectorySelection/Concordance.lean`, `statistical_significance` | Returns constant `0.05`, independent of the result. | Derive or implement a justified test with null assumptions and calibration. |

These definitions are outside the proof ledger's syntactic inventory but remain explicitly provisional. Downstream work must not cite them as evidence of replication, significance, or power.

## Audit recipe

From the repository root, search Lean source for `sorry`, `axiom`, and `native_decide`, including comments and declarations, then inspect each match against this ledger; exclude library dependencies under `.lake/`. Compile the full library with `lake build` (or `cd lean && lake build`) using the pinned toolchain/lockfile. A passing build alone does not establish that the ledger is complete or that the science is sound.
