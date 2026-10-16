# ARG phase: progress and plan (2026-10-15)

Companion to `docs/CURRENT-STATE.md` (ADDENDUM 2026-10-15). Short form: what exists,
what the next block of work is, and how each item is judged.

## Why the ARG (one paragraph)

Window/cluster methods could not resolve local ancestry at any array density
(run-length ratio 1.00-1.05 across densities for moderns and ancients; see
CURRENT-STATE). The ARG removes windows entirely: ancestry is the topology, time is
tsdate, and every locus statistic becomes a tskit query. Measured cost: ~10 min per
chromosome at ~1,300 samples; whole genome in parallel ~15 min; ~4 GB on disk.

## Status of components

| # | component | state | evidence |
|---|---|---|---|
| 1 | ARG inference on raw (unimputed) AADR calls, radiocarbon-dated | **done** | `analysis/arg_build.py`; 41k trees, chr1 |
| 2 | Modern clock anchors (fix oldest-tip age pinning) | **done** | mixed cohort 1,000 ancients + 322 EUR moderns |
| 3 | Dense merged substrate: AADR 1240K ancients + Akbari 1000G moderns, position-matched | **done** | 1,077,154 shared sites; `analysis/arg_merge.py` |
| 4 | Allele **origin-age** test vs frequency | **done, and superseded** | detects nothing at LCT: z=-0.63. Reason: allele age != selection time (LCT allele is old; the *rise* is recent) |
| 5 | **Frequency-trajectory** estimator (CLUES-style) | **this block** | subject of the work below |
| 6 | Trajectory validation on positive controls (LCT, SLC45A2, SLC24A5, FADS) | **blocked on 5** | must show a recent rise although the alleles are old; if not, 5 is wrong |
| 7 | AMY1 read-off on the validated scale | **blocked on 6** | no prediction; report position in the control-calibrated distribution |
| 8 | Punctuated test: is the rise sharp, single-background, era-dated | **blocked on 6** | the project's motivating question |
| 9 | Genome-wide scan (all 22 chromosomes, merged substrate) | **blocked on 6** | ~15 min once the statistic is validated |

## Item 5/6 design (registered before running)

Estimator: for a variant with derived-allele carriers, compute the **carrier-lineage
fraction** on the dated ARG:

    rho(T) = (# edges present at time T whose child lies in the carrier subtree)
             / (# edges present at time T)

`rho` is the proportion of ancestral lineages carrying the derived allele at time T —
the standard ARG frequency proxy (Relate-style). Reported per variant:

- present-day frequency (from genotypes),
- `T_rise` = oldest time at which rho >= 50% of its present value (when the allele
  started rising, as opposed to when it first existed),
- `width` = T(25%) - T(75%): small width = punctuated rise, large = gradual,
- full rho(T) curve for inspection,
- the same for the frequency-matched window null.

Controls (must pass before AMY is read): LCT rs4988235 (chr2:136,608,646), SLC45A2
rs16891982 (chr5:33,951,693), SLC24A5 rs1426654 (chr15:48,426,484), FADS (chr11:61.6 Mb
block). Prediction: T_rise in the last few thousand years (Holocene), i.e. far more
recent than the allele's origin, and a narrow width.

Falsification: if the controls show no recent rise, the estimator (not the biology) is
the problem and item 5 is revised before anything else.

---

## RESULTS (2026-10-15): item 5 implemented, item 6 partially validated

### What was built

- `analysis/arg_traj.py` — carrier-lineage-fraction trajectory on the **local tree** at a
  site: rho(T) = (edges present at time T whose child is in the carrier subtree) / (edges
  present at time T). Reports present frequency, T_rise (oldest time with rho >= 50% of its
  present value), width, and the curve.
- `analysis/arg_modern.py` — **modern-only ARG** (503 EUR 1000G, all time-0 samples):
  contemporaneous sampling makes the lineage-fraction trajectory interpretable (no
  temporal-sampling confound) and enables tsdate's inside-outside method. Unphased
  heterozygotes are resolved at random (standard workaround; seeding fixed).
- `analysis/arg_validate.py` — per-window report of carrier-clade age for target variants
  and their rank among common variants (freq >= 0.30).

### Key result: the estimator lights up at LCT

Modern ARG, chr2:135.5-137.5 Mb, 503 EUR:

| variant | freq | carrier-clade age | rank among common (n=172) |
|---|---|---|---|
| **LCT rs4988235 (chr2:136,608,646)** | 0.505 | **35.6 gen (~1,033 yr)** | **4** |
| window common-variant median | — | 516 gen | — |

And the window shows a **cluster** of young common variants at 136.83-136.92 Mb (ages
30.8-35.6 gen) — the linked haplotype block a recent sweep produces. Contrast: an
**ancestors-only** ARG puts rs4988235 at z = -0.63 (no signal), and the older
allele-**origin** statistic is blind to it by construction. So: punctuated rise detected,
gradualist/origin statistics blind — exactly the distinction this project turns on.

### The limitation that now binds (next action, item 9')

The merged substrate gives only ~300-350 sites per 2 Mb window after QC (1240K x Akbari
intersection), which is too sparse for fine clade resolution: SLC45A2 rs16891982 (freq 0.93)
comes out old (1,260 gen, rank 130/200), and per-window "youngest common variant" ages are
similar across AMY (39-70 gen), FADS (61-81), SLC45A2 (35-48) and LCT (31-36), i.e. the test
is not yet discriminative.

**Next action**: rebuild the modern ARG on the **full dense variant set** (the Akbari
AllVariants sites) rather than the 1240K intersection. For moderns those genotypes are real
1000G WGS, so a dense modern ARG is both valid and 10-50x better resolved (thousands of
sites per window), while the ancients stay on real AADR 1240K calls as an independent
frequency check across time. Then re-run LCT / SLC45A2 / SLC24A5 / FADS controls and only
then read AMY1.

---

## Assessment (2026-10-15): what this metric gives, and what it does not

### Provenance — yes, borrowed (state it plainly)

Allele age from a coalescent tree is classic (Slatkin & Rannala 2000; Griffiths &
Tavare 1998; Thomson 2000). ARG-scale versions: Relate allele age and omega (Speidel
et al. 2019), Albers & McVean 2020, Fan et al. 2019. Frequency-trajectory / selection
from a dated ARG: CLUES (Stern, Wilton & Nielsen 2019). Our `arg_traj.py` is a
simplified **carrier-lineage-fraction proxy** — essentially Relate's frequency curve
without the probabilistic model. It is not a new method, and using established
estimators is the right call: results stay comparable to the literature.

### What is genuinely ours

1. **Real ancient genotypes as an independent temporal check.** Most ARG-selection work
   has moderns only, so the inferred trajectory cannot be tested. We hold raw,
   unimputed AADR calls for ~1,000 dated individuals and can ask whether the ancient
   allele frequencies at their own dates follow the trajectory inferred from moderns.
2. **ARG-derived backgrounds for the carriage/selection question.** Shared-branch-length
   relatedness on a dated ARG is a real ancestry partition with a time axis, unlike the
   ELAI per-window clusters (which failed the tract test at every density).
3. **The punctuated framing made testable**: a rise is a sharp event at a datable time
   on a specific background — the quantity the gradualist GLM absorbs.

### What it gives us so far (concretely)

- Working, fast pipeline: dated ARG on real ancient data, genome-wide in ~15 min.
- LCT rs4988235 detected as common-but-young (35.6 gen, rank 4/172) with a linked
  young haplotype cluster; the ancestors-only and allele-origin analyses are blind to it.
- A calibrated limitation: ~300 sites/window is not discriminative.

### Honest limitations of the metric

1. **Proxy, not model**: no uncertainty quantification, no sampling-time correction
   (dodged by using moderns only), no recombination-aware prior. CLUES exists; we should
   port/use it or label ours explicitly as a screen.
2. **Crude null**: rank among common variants in the same window ignores local mutation
   and recombination rate, gene density, background selection, and multiple testing.
3. **Ages may be systematically young** (LCT ~1,000 yr vs classical ~5,000-10,000 yr
   estimates): needs calibration against variants with independently known dates.
4. **Sparsity quantizes clade ages** (many linked variants share one node age), which
   also inflates apparent clusters.
5. **Not yet tied to Akbari's X** in a principled way; the comparison is still ad hoc.

### TODO (next, in order)

1. **Dense modern ARG** — rebuild from the full AllVariants sites (real 1000G WGS, 10-50x
   more sites/window); removes the resolution limit. Keep ancients on real AADR calls.
2. **Proper null** — frequency-matched AND locally mutation/recombination-matched;
   genome-wide, with multiple-testing control.
3. **Age calibration** — check the estimator against variants with independently known
   ages before believing any absolute date.
4. **Controls** — LCT (passing), SLC24A5, FADS, SLC45A2 must all behave sensibly.
5. **AMY1** — read only after 1-4.
6. **Carriage/selection decomposition** — ARG-background clustering + trajectory timing.
7. **Repo hygiene**: this repo has NO git remote (117 commits are local-only); the ELAI
   patch is the only public artifact. Add a remote and push if the work should be visible.
