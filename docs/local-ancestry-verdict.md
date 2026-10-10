# Local-ancestry verdict: the states are not identified

> Parent-run record. Closes the local-ancestry route with numbers. Upstream artefacts:
> `docs/elai-amy-long-run-results.md`, `docs/elai-settling-sweep.md`,
> `docs/local-ancestral-state-method.md`.

## What was run

Ten-arm settling sweep (4,013 s wall, 10 arms parallel) followed by three full-density fits
with the settled configuration: **K=4 at 100 steps, seeds 20261013 and 20261014**, plus
**K=5 at 100 steps, seed 20261015** (still running at the time of writing). Panel:
54,370 retained markers, N=1,000 dated ancients, AMY core excluded, virtual site at
chr1:104,121,000.

## Finding 1: the convergence verdict is decided by numerical noise

| fit | steps | status | first qualifying | final relative increment | 2nd-half mean abs | negative increments |
|---|---:|---|---:|---:|---:|---:|
| K=4 seed 20261014 | 100 | **PASS** | 48 | +1.20e-05 | 4.65e-05 | **0/50** |
| K=4 seed 20261013 | 100 | **FAIL** | 56 | -1.61e-05 | 5.25e-05 | **12/50** |

Two seeds, identical configuration, opposite verdicts. The increments have decayed to the
**noise floor** (~5e-05 relative on a likelihood of ~2.5e6, i.e. ~50-100 absolute in 2.5
million), which sits only ~2x below the 1e-04 tolerance. So the convergence criterion is
measuring floating-point noise rather than a clean decay, and the earlier 200-step runs
failed precisely because they ran past the noise floor and the rule requires the *final*
window to qualify.

[inference] This is a **criterion-validity** problem, not merely a fit problem: at this
tolerance, the verdict is set by noise.

## Finding 2 (the decisive one): the site posterior is unidentified and irreproducible

The convergence criterion is only a proxy. The quantity the analysis actually needs is the
local haplotype-sharing state at the excluded AMY core. Comparing the two K=4 seeds at
chr1:104,121,000 over all 1,000 individuals x 16 ordered pairs:

| diagnostic | value | reference |
|---|---|---|
| ordered posterior sums per individual | **exactly 1.00000** (both seeds) | correct: the patch works |
| mean abs difference between seeds | **0.079** | values average 0.0625, so the difference is as large as the values |
| correlation between seeds | **-0.19** | not positive |
| dominant state-pair agreement | **3.4%** | chance is 6.25% (1 of 16) - **worse than chance** |
| dominant unordered-pair agreement | 4.5% | chance 10/16 |
| median per-individual max posterior | 0.18 (A) / 0.28 (B) | uniform is 0.0625 |
| individuals with max posterior > 0.5 | **0.4% (A) / 5.1% (B)** | a confident assignment would be common |
| median posterior entropy | 2.42 (A) / 2.14 (B) | uniform ln16 = 2.77 - i.e. nearly flat |

**[inference] The local sharing states are not identified.** The posterior is close to
uniform (entropy 2.1-2.4 against a maximum of 2.77), almost no individual receives a
confident state (0.4-5.1% above 0.5), and two seeds agree on the dominant state pair
*less often than chance would*. Whatever the likelihood is doing, the state assignment we
would need to condition on is not being determined by the data.

This also explains why the two seeds disagree on convergence: they are settling into
different regions of a flat ridge, not converging to a common optimum.

## Finding 3: the knobs available do not fix it

The sweep already swept what the model exposes: `S` (1/2/3), `K` (2/3/4/5), mixture
generations `-mg` (5/10/100/1000), warm-up `-w` (0/1), and step count (50-200). Results:
`-mg` is **inert** (identical statistics at 5/10/100/1000, so the admixture-generation
parameter does not enter this fit), `-w 1` makes settling worse (25.7% wrong-sign), step
count controls the verdict, and `K=3/4` settle more readily than `K=5` - but none of that
produces an identified state assignment.

## Verdict

**The local-ancestry conditioning route is closed on this data and configuration.** It was
not closed by the reference-panel coverage problem (attempt 1), nor by the wrong model class
(attempt 2), nor by export size - those were all fixed, and the patched tool runs at full
density with correct arithmetic. It is closed because **ELAI does not resolve the local
haplotype-sharing states on this panel**: the posterior is flat and seed-unstable.

This is consistent with the published caveats the design lane found: local-ancestry
agreement between methods is only 0.08-0.50, and Mies & Mathieson (2026) note their own
frequency deviations "support either within-background change or underestimated ancestry
shifts".

## What would be needed to reopen it

A different pipeline, not a tuning of this one: principled **phasing** of the ancient
genomes followed by **RFMix or Gnomix** on **curated, adequately sized reference panels** -
which is what the published Neolithic study used, and whose own inter-method local agreement
was still only 0.08-0.50. That is a separate project, and it would have to clear the same
structure and stability gates before any locus could be conditioned on it.

## Consequence for the project

The temporal question - does the AMY allele change **within** a fixed local ancestral state
- therefore remains **untested**. The project's established results are unaffected and stand
on their own: AMY sits in the discarded band (`|X|` = 3.795 against the 5.45 cut), the
binarised short/long allele clears FDR-50% (X-equiv 3.64) and survives the structure
correction (~97% retention), the signal is ~94% carriage on the global-label decomposition,
and there is no low-`|X|`/high-replication class genome-wide.

## Confirmation at K=5 (100 steps, seed 20261015)

| diagnostic | K=4 seed A | K=4 seed B | K=5 seed 20261015 |
|---|---|---|---|
| convergence status | FAIL | PASS | **FAIL** (3/50 negative = 6%, limit 5%) |
| final relative increment | -1.61e-05 | +1.20e-05 | +1.21e-05 |
| ordered posterior sum | 1.00000 | 1.00000 | 1.00000 |
| median per-individual max posterior | 0.181 | 0.281 | 0.200 |
| individuals with max > 0.5 | 0.4% | 5.1% | **2.3%** |
| median entropy | 2.42 | 2.14 | 2.54 |
| entropy / uniform | 0.87 | 0.77 | **0.79** |

**[inference] More clusters do not help.** K=5 fails convergence by 6% versus a 5% limit -
again decided at the noise level - and its posterior is just as diffuse relative to its
larger state space (entropy ratio 0.79 versus K=4's 0.77-0.87), with 97.7% of individuals
lacking any state above 0.5. The states remain unidentified at both K=4 and K=5.

**The verdict stands, and is now supported at two cluster counts:** the local
haplotype-sharing states are not resolved by ELAI on this panel, so the within-state
temporal test cannot be built on them.
