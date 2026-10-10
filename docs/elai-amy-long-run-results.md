# ELAI dense long-run results (gated STOP)

> Parent-run record. Three parallel fits at full density were launched by the parent and
> this file records their outcome, the gate verdict, and the diagnosis. The run outputs
> live under `data/derived/elai_amy_dense/work_long/` (gitignored); the numbers below are
> the durable record.

## What ran

| fit | config | wall clock | steps | complete |
|---|---|---:|---:|---|
| K=5 seed A | `-p 1 -C 1 -c 5 -w 0 -s 200 -mg 10 -R 20261010 --pkk-site 1:104121000 --pk1` | 11,845 s (3.29 h) | 200 | yes |
| K=5 seed B | same, `-R 20261011` | 11,931 s (3.31 h) | 200 | yes |
| K=3 seed A | same, `-c 3 -R 20261010` | 5,543 s (1.54 h) | 200 | yes |

Input (persistent, built once): **54,370 retained markers**, fitted **N=1,000** dated
ancients, **9,854** independent PC markers. Marginal export: **271,850,000** values per
K=5 fit, drained through a FIFO, never persisted.

## Gate verdict: FAIL (no biological use authorized)

| gate | status | detail |
|---|---|---|
| exclusion | **PASS** | 0 retained core markers, 163,110 checked; site 104,121,000 inside core 103,915,000-104,327,000 |
| virtual_site.arithmetic | **PASS** | ordered-sum error 4.4e-16; q symmetry 0.0; triangle error 6.7e-16 |
| coverage.usable_input | **PASS** | 1,000/1,000 usable |
| coverage.marginal_export | **PASS** | complete, 271,850,000 values, norm error 1.0e-3 |
| em_convergence | **FAIL** | see below |
| known_structure | **NOT TESTED** | reason: no converged primary, labels not accepted |
| stability | PARTIAL | seed / cross-K / virtual-argmax / thinning all NOT TESTED |
| coverage (joint) | PARTIAL | full ordered-joint coverage NOT TESTED (export-size driven) |
| virtual_site.integration | NOT TESTED | frozen-linkage replays not run |

`verdict = FAIL`, `fits = 3`, `biological_use_authorized = false`.

## Diagnosis: the K=5 EM does not settle

Convergence trace statistics (tolerance 1e-4 relative):

| fit | first qualifying step | second-half mean abs increment | second-half negative increments |
|---|---:|---:|---:|
| K=3 seed A | 50 | **8.31e-06** | 13/100 |
| K=5 seed A | 46 | **2.00e-05** | 41/100 |
| K=5 seed B | 50 | **2.15e-05** | 33/100 |

**The oscillation is reproducible across seeds.** Both K=5 fits reach the 1e-4 tolerance
by step ~46-50 and then oscillate for the remaining 150 steps at ~2e-5 relative with
30-40% of increments going in the *wrong* direction. The likelihood does not decay toward
a fixed point; it bounces on a plateau. Representative K=5 seed A tail (absolute
increments on a likelihood of ~2.54e6): `+101, -50, +3, -41, +31, -14`.

So the gate verdict is **substantive, not an artifact**: a sign-flipping increment at
constant magnitude is not convergence, and withholding labels was correct. (For K=3 the
oscillation is smaller, ~8e-6 with 13% negatives, which is closer to settling but still
fails the preregistered 3-consecutive-non-negative criterion.)

An earlier parent note describing this as a gate artifact is **retracted for K=5**. The
non-negativity condition detects real non-settling here; it is not merely floating-point
jitter.

## Likely cause to test next

The configuration is `S=1` (`-C 1`), i.e. a single upper state, chosen as the linked-model
limit so that the lower `K` haplotype-copying clusters carry the sharing state. With a
single upper state the upper-level parameters have little to identify them, which is a
plausible source of a flat likelihood ridge and hence oscillating EM. Candidates to probe,
all cheaply, in a short EM-settling sweep (30 steps, several settings in parallel):

- upper states `S` = 1 / 2 / 3;
- lower clusters `K` = 3 / 4 / 5 (K=3 settles markedly better than K=5);
- prior scale `-mg` = 5 / 10 / 20;
- warm-up `-w` = 0 / 1.

The decision rule for that sweep is whether the second-half mean abs increment **decays**
with step count and whether the negative fraction falls below ~5%, i.e. whether the fit
settles rather than merely dipping under the tolerance.

## What this run established, and what it did not

**Established:** the patched binary runs on real imputed data at full density; the
exclusion and virtual-site arithmetic are exact to machine precision; the marginal export
path works at scale; and **ELAI's K=5 EM does not converge on this data/configuration**,
reproducibly.

**Not established:** nothing about ancestry structure, cluster stability, or any allele
frequency at the AMY locus. Because labels were withheld, no local-state summary exists,
and the temporal analysis remains unstarted.

This is the honest ceiling of a marginal-export fit: even with convergence, the missing
mandatory checks (full ordered-joint coverage, frozen-linkage integration) yield at best
PARTIAL — and here convergence itself failed, so the verdict is FAIL.
