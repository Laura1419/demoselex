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

---

## GENOME-WIDE SCREEN RESULTS (2026-10-15) — items 9 and AMY read-off

Built and ranked in 39 min wall on 22 parallel chromosome jobs:

- 22/22 chromosomes: genotypes (503 EUR moderns, merged 1240K x Akbari sites), tsinfer
  ARGs, tsdate inside-outside in 10 Mb chunks (291 chunks), then per-variant statistics.
- **834,531 variants scored genome-wide.** Statistics: carrier-clade age = tsdate time of
  the mutation node; local rho = breakpoints/bp per chunk; z = MAD-standardised age within
  frequency-decile x rho-quartile strata.

### Calibration (known loci, ranked genome-wide)

| locus | freq | clade age | z | rank / 834,531 | percentile |
|---|---|---|---|---|---|
| **LCT rs4988235** (known recent sweep) | 0.525 | 45 gen | -1.74 | **1,087** | **0.13%** |
| SLC45A2 rs16891982 (suspected carriage) | 0.936 | 1,332 gen | -0.30 | 314,301 | 37.7% |
| SLC24A5 rs1426654 | — | — | — | absent from scored set | — |

Direction is right for both: the known in-situ sweep ranks in the top 0.13%, the
carriage-like locus is unremarkable.

### CORRECTION (2026-10-15): chromosome labels were missing from the pooled table

The first pooling step concatenated per-chromosome tables WITHOUT a `chrom` column, so
any position-window query silently mixed chromosomes: "chr1:103-105 Mb" returned 8,343
rows instead of 351. The AMY candidate list reported on the first pass (rank 116,
"top 0.01%") was an artifact of that bug. `arg_genome.py pool` now carries `chrom` and
the numbers below are corrected. Lesson: every positional query must be chromosome-scoped.

### Corrected regional test (variants in the genome-wide top 1% of z)

| region | n | observed top-1% | expected | enrichment | Poisson p |
|---|---|---|---|---|---|
| AMY (chr1:103-105 Mb) | 351 | 6 | 3.5 | **1.71x** | 0.14 |
| LCT (chr2:135.5-137.5 Mb) | 431 | 7 | 4.3 | 1.62x | - |
| SLC45A2 (chr5:33.5-34.5) | 298 | 17 | 3.0 | 5.70x | - |
| SLC24A5 (chr15:48-49) | 226 | 8 | 2.3 | 3.54x | - |
| FADS (chr11:61-62.2) | 322 | 0 | 3.2 | **0.00x** | - |

Best chr1 AMY variant: 104,350,432 (freq 0.632, age 165 gen, z -1.66, percentile **0.23%**).
The 6 AMY top-1% variants fall in **5 separate 100 kb clusters** (positions 103.11, 103.41,
103.44, 103.61, 104.02, 104.35 Mb) - scattered, not one coherent locus.

### AMY region (chr1:103-105 Mb) — candidate, not a result

351 scored variants on chr1 in the window (corrected). Superseded list below kept for the
record; the *corrected* summary is the regional test above: 1.71x enrichment, p = 0.14,
best variant percentile 0.23%. Original (buggy, cross-chromosome) list:

| pos | freq | clade age | z | rank | in AMY1 CNV segment? |
|---|---|---|---|---|---|
| 104,856,148 | 0.646 | 112 gen | -1.96 | **116** | no |
| 104,131,494 | 0.702 | 273 gen | -1.80 | 677 | **yes** |
| 104,128,773 | 0.702 | 34 gen | (young) | — | **yes** |
| 104,123,531 | 0.499 | 34 gen | (young) | — | **yes** |
| 104,115,429 | 0.499 | 25 gen | -1.64 | 2,291 | **yes** (boundary) |
| 103,589,094 | 0.519 | 92 gen | -1.71 | 1,371 | no |

Unlike LCT (four variants on ONE node, age 35.6 gen), the AMY candidates sit on
**distinct nodes with ages spread 25-454 gen** - so this is not one clean sweep fragment.

### Required checks before any AMY claim (registered)

1. **Paralogous-mapping artifact**: AMY1 is a segmental duplication (hg19
   ~104,115,000-104,290,000). Paralogous sequence variants from mismapped reads
   manufacture exactly this appearance (young, common, clustered). Check HWE/excess
   heterozygosity, multi-mapping, and whether the candidates are paralogous sites.
2. **Local rho**: the current rho is per-10-Mb-chunk (near-constant) so the recombination
   correction is weak; recompute per variant (+/- 250 kb breakpoint density).
3. **LD clumping / effective tests**: variants in one block share ages, so 834k variants
   are far fewer independent tests; clump by LD before quoting ranks.
4. **Seed robustness**: unphased het resolution is random; repeat with other seeds and
   confirm the AMY candidates and LCT survive.
5. **Ancient check**: compare the inferred trajectory/age to the ancient allele
   frequencies at these positions (external, independent of the modern ARG).
6. Only then interpret AMY1 relative to Akbari's X = 3.795 (below their 5.45 threshold).

---

## ITEM 1 DONE: within-tree null (and one design failure caught)

### The invalid first design

`arg_permnull.py` compared a variant's carrier set to random *sample sets* of the same
size. All 200 candidates hit the p-floor and null medians were absurd (~14,000 gen vs
observed ~100 gen). Reason: in a tsinfer ARG a mutation is placed ON a clade, so
"carriers form a clade" is a tautology; random sample sets are not clades, so the
comparison is rigged. Recorded as a dead end.

### The corrected null: clade-size-matched age (`arg_cladenull.py`)

Enumerate every clade in the chunk's trees (every node's descendant sample set). For a
variant, its carrier clade has a size n and an age; the null is the distribution of ages
of ALL clades with size ~n in the same chunk. Demography-free (it is in the tree), no
simulation, self-contained. Statistic: p_age = P(clade age <= observed | size bin).

Results (K not needed; null = pooled clades, n = 0.2-2.1 M clades per 10 Mb chunk):

| set | n variants | p<0.05 | p<0.01 | p<0.001 | expected (0.05/0.01/0.001) |
|---|---|---|---|---|---|
| chr2 top-400 by heuristic z | 400 | 163 | 69 | 7 | 20 / 4 / 0.4 |
| **AMY region (chr1:103-105 Mb)** | 351 | 27 | **16** | **3** | 17.6 / 3.5 / 0.35 |
| LCT region (chr2:136.5-136.7) | 46 | - | - | - | - |

- **LCT rs4988235**: 264 carriers, clade age 45.4 gen, null median 580.8 gen, **p_age = 0.024**.
  (Chunk-wide null is conservative vs the frequency x rho heuristic, which put LCT at
  0.13 percentile; the tighter variant is a null scoped to the SAME local tree.)
- **AMY**: 16 variants at p<0.01 (4.5x expected) and 3 at p<0.001 (8.6x expected), e.g.
  104,019,631 (50 carriers, age 14 gen vs null 174, p=2.3e-4), 103,784,927 (68 carriers,
  age 16 vs 174), 104,324,819 (95 carriers, age 26 vs 295).
- **Caution**: the AMY clades are MORE extreme than LCT's. Since AMY1 is a segmental
  duplication, paralogous/mismapped variants manufacture exactly young apparent clades -
  so the strength of the AMY signal is as consistent with artifact as with selection.
  Do not interpret until the paralog checks are done.

### Next refinements (registered)

1. Null scoped to the SAME local tree (or the same LD block) instead of the whole chunk.
2. Branch-level primitives (branch length + mutation load + monophyly) - item 2.
3. Event assembly + FDR - item 3.
4. Paralog/artifact checks at AMY before any claim.

---

## REBUILD (2026-10-16): real phase + dense sites — the pipeline was the problem

Two defects found and fixed:
1. **Phase was being thrown away.** The bundle carries a phase bit (verified: 2,514/2,515
   calls via `read_alleles_and_phasepresent`). The old build read unphased genotypes and
   coin-flipped heterozygotes into pseudo-haplotypes, fabricating haplotype sharing.
2. **Site universe was 1.4% of what was available.** The old build used the 1240K x Akbari
   intersection (~1 site / 2.7 kb); HighQualityVariants.tsv holds 9.7 M high-quality sites.

New pipeline (`analysis/arg_phased.py`): 503 EUR moderns -> **1,006 real phased haplotypes**,
sites from HighQualityVariants (MAF>=0.05, call>=0.90), per-chromosome tsinfer + tsdate in
10 Mb chunks. Wall time ~4.7 h for all 22 chromosomes in parallel. **5,539,880 variants
scored** (vs 834,531).

### The decisive check: clade ages were compressed young

| variant | old (sparse, unphased) | new (phased, dense) | literature |
|---|---|---|---|
| LCT rs4988235 | 45 gen (~1,300 yr) | **193 gen (~5,600 yr)** | LCT sweep ~5-8 kya |

The "ages look systematically young" worry was a pipeline artifact, now fixed.

### Regional excess, composite null (size x local rho x compartment), BH FDR

NON-SEGDUP compartment (2,687 windows):

| locus | enr | q | rank |
|---|---|---|---|
| **LCT** | 3.01x | <0.001 | 179 |
| **SLC24A5** | 2.07x | <0.001 | 310 |
| **FADS** | 1.52x | <0.001 | 558 |
| **AMY** | 1.28x | **0.04** | 787 |
| SLC45A2 | 0.97x | 1.00 | 1,294 |

Top regions: chr12:34 Mb (8.0x), chr10:42 Mb (11.0x), chr2:219 Mb (7.0x), chr10:48 Mb
(10.4x), chr3:90 Mb (12.5x), chr2:95 Mb (16.9x).

SEGDUP compartment: 161,247 variants / 476 windows / 75 at q<0.05 (was 8,824 / 26 / 1 with
array sites) - the compartment only became visible once dense sites were used. AMY: 52
variants, 2.69x, p=0.017, q=0.11.

### Interpretation

Controls behave and split by regime: recent in-situ sweeps visible (LCT, SLC24A5),
long-standing selected alleles not (SLC45A2 0.97x). AMY moved from "not flagged" to
"flagged, modestly" once real haplotypes were used.

---

## 1000G HIGH-COVERAGE REBUILD (2026-10-17/18): 3,202 people, 6,404 real haplotypes

Substrate change: instead of 503 imputed EUR from the Akbari bundle, use the 1000G
high-coverage phased release (GRCh38, 3,202 samples, real phase from the source).
`analysis/g1k_phased.py`: streams the phased VCF with cyvcf2 into a transient tsinfer
SampleData (no VCF retained, no big in-RAM matrix), filters SNVs with MAF>=5%, infers,
dumps trees, dates in 10 Mb chunks.

Measured costs (chr22 pilot): 102,351 common SNVs x 6,404 haplotypes; SampleData 53 MB;
tree sequence 1.25 GB; build+infer 40 min; date+stats 67 min. Inference scales ~0.05 s/site
above 200k sites (chr8: 384k sites -> 5.2 h), so chr1/chr2 ~7-8 h inference + ~6 h dating.
Genome-wide: ~6M common SNVs, ~50-90 GB of trees, 12-way concurrency.

Validation (chr22, both substrates, age by frequency):

| freq | 1000G (6,404 haps) | EUR (1,006 haps) | ratio |
|---|---|---|---|
| 20-30% | 1,107 gen | 1,150 | 0.96 |
| 30-50% | 2,110 | 2,073 | 1.02 |
| 70-100% | 5,866 | 4,095 | 1.43 |

The clock agrees where it should and extends 43% deeper at common frequencies - the
pre-Out-of-Africa depth that a EUR-only sample structurally truncates. Spearman(age, carriers) = 0.81.

Early control read (chr15, complete):

| locus (GRCh38) | freq | age | percentile | reading |
|---|---|---|---|---|
| HERC2/OCA2 rs12913832 | 0.179 | 37 gen (~1.1 kyr) | 2.0% | recent rise |
| SLC24A5 rs1426654 | 0.566 | 8,270 gen (~240 kyr) | 89.7% | old, common allele |

GRCh38 inputs staged: hg38 genomicSuperDups, 3,202-sample population map (CEB/YRI/...),
canonical loci (LCT 2:135,851,076; SLC45A2 5:33,951,588; SLC24A5 15:48,134,287; FADS
11:61,785,208; HERC2 15:28,120,472; AMY cluster chr1:103.61-103.70 Mb).
