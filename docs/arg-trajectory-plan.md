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
