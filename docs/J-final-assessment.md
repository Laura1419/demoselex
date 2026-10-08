# J. Final assessment

> Capstone for the project. Verdict first, then the evidence chain, then what would
> overturn it. Labels: **[observation]** measured, **[inference]** conclusion drawn,
> **[speculation]** unverified.

## Verdict

**[inference] The initial criticism is NOT SUPPORTED.** On the decisive like-for-like
data (Akbari's own individuals and variant set), there is no detectable class of loci
whose signal is weakened by the population-structure correction while showing replicated
directional change across population histories. The hypothesised low-Akbari /
high-replication class is empty at the power we can reach.

This is stated with a specific scope, not as a claim that Akbari's method is perfect:

- What is **falsified** is the specific prediction that replicated within-history or
  between-history directional change identifies a set of loci that Akbari's correction
  suppresses.
- What is **established** is that the two statistics *agree*: loci with a strong Akbari
  signal are ~8.6x more likely to also show high within-trajectory replication.
- What is **limited** is power: with `N_eff ~ 4-6` independent histories, several tests
  cannot reach genome-scale significance by construction, so "no evidence" is not the
  same as "no such loci anywhere".

## Evidence chain

### 1. The three regimes were tested, one at a time, with pre-registration where it mattered

| regime | statistic | data | result |
|---|---|---|---|
| within-history change (pooled across histories) | `trajectory_statistic`: per-history slopes, concordance, relatedness-discounted GLS metal-analysis | Tier-1 AADR (5,788 we individuals, 1.23M sites) | low-\|X\| pass rate **0.0427** vs 0.05 null; high-\|X\| **0.0924** (~2.1x) |
| **same, decisive data** | same | **Tier-2: 5,737 individuals, paper's imputed 9.7M-variant QC set** | **Q3 = 0**; low-\|X\| **0.0438** vs 0.05 null; high-\|X\| **0.3763** (**8.6x**) |
| selection in a *subset* of histories | pre-registered max-over-subsets (max \|z\| over subsets of size 2-4) | Tier-1 panel, 2,000 candidates | **0 pass BH-FDR** (min q 0.639); enrichment 1.13x |
| between-history ordered rise (histories internally flat, successive turnovers higher) | pre-registered relatedness-aware GLS of history mean frequency on history placement | Tier-2 panel, 8,000 variants | null **calibrated** (control pass 0.044); low-\|X\| **0.044** = null rate; high-\|X\| 0.087; **0 pass BH-FDR**; prediction **falsified** |

All three pre-registered predictions point the same way: enrichment is at the
Akbari-significant end, not the low-Akbari end.

### 2. Simulation says where the method *should* win, and it is not where the data are

| selected histories | Akbari GLMM power | max-over-subsets power |
|---|---:|---:|
| 1 of 6 | 0.25 | **0.95** |
| 2 of 6 | 0.80 | **1.00** |
| 6 of 6 | 1.00 | 1.00 |

With two measured costs: calibrated FPR under pure replacement **0.112 vs 0.045**, and
the uncalibrated max device runs at **0.55-0.64** FPR. So the advantage exists, is
narrow, and is confined to *heterogeneous* selection - and that regime is not evidenced
in West Eurasian data at this power.

### 3. The structural reason the hypothesis fails

**[inference]** The two methods' blind spots **overlap**:

- Akbari is blind to the component of the time vector collinear with the dominant
  ancestry axis (formalised: retained information = `(1 - cos^2 theta) + cos^2 theta *
  sigma^2/(sigma^2 + V_mode)`; the slope is unidentifiable at full alignment).
- Our statistic is blind to *carried-in* change: when an expanding population brings the
  allele, every history is internally flat, so within-history slopes are ~0.

But "collinear with ancestry" and "carried in" are the same physical situation. So in
precisely the regime where Akbari loses the signal, our statistic has nothing to measure.
The pooled GLMM is also the *efficient* estimator for homogeneous change, so where the
signal is clean (all histories rising) Akbari wins outright.

**[inference]** The honest reframing that came out of this work: the method is a
**heterogeneity detector**, not a better trend detector. What costs Akbari is not the
structure correction per se - that correction is doing its job - it is the single pooled
slope.

### 4. What the formalization contributed

- `Replacement.lean`: under a replacement-only model the concordance statistic has
  **exactly** the exchangeable null distribution, and `P(perfect pair concordance) = 1/2`.
  So rejecting that null is meaningful, and replacement alone cannot manufacture evidence.
- `Statistic.lean`: `N_eff` is **duplication-invariant** - adding a copy of an
  already-represented history does not increase the effective count.
- `Collinearity.lean`: the GLMM's loss of identifiable information is a theorem, not a
  hand-wave, with an explicit formula.
- **No `sorry`s and no project axioms**; `lake build` green (3103 jobs). Caveat:
  `lean/AXIOMS.md` still lists three inherited placeholder definitions
  (`concordance_score` returns a constant `0.5`, `check_replication` uses a `0.7` cutoff,
  `statistical_significance` returns a constant `0.05`) which must not be cited.

### 5. Data and reproducibility

- Akbari per-variant statistics (674 MB) and the processed imputed PLOS2/PLINK2 bundle
  (26.3 GB, md5 verified), AADR v66.p1 1240K + HO (11 GB), AMY tables - all downloaded
  with recorded sha256/md5, logged in `docs/artifacts/download_log_2026-10-08.tsv` and
  `docs/tier2-data.md`.
- Every run has a provenance manifest; the task graph, decision log, and journals are in
  version control (90+ commits).
- Known data pitfalls recorded: `X` is `Z/sqrt(CF)` with `CF = 2.28` and is *not*
  centred (mean +0.412), so "low Akbari" must be defined by absolute thresholds
  (`|X|` 3.16 = FDR 50%, 5.45 = pi>99%), never by `X ~ 0`; the imputed bundle stores
  integer-like hard calls; the interrupted Tier-2 resume appended past the object end and
  required truncation to the API-listed size.

## Limits - read before quoting the verdict

1. **Replication power ceiling.** `N_eff ~ 4-6`. The between-history ordering test has a
   `~1/m!` floor before multiple testing; genome-scale BH-FDR is **unreachable by
   construction** with the finite matched null used here (P floor 0.0024 vs required
   1.7e-5 over 3,000 tests). "0 pass FDR" is therefore largely a resolution statement.
2. **The replication units are imperfect.** The 6 merged groups contain a
   10,058-individual `G0` bloc that is an over-merged period/region artefact; the honest
   distinct-history count is ~4, and there is no natural break in the relatedness
   spectrum.
3. **Tier-2 coverage.** 5,737 of 11,934 trajectory individuals are in the imputed bundle;
   12,000-8,000 of 9.7M variants were scored. Hard-call dosages mean imputation
   uncertainty is not propagated.
4. **The AMY case study is incomplete.** 537 ancient carriers, 522 with dates and
   coordinates, 288 across 10 structural haplotypes - but the five case-study questions
   are blocked on the trajectory linkage and sequence-level data. The mutation-level
   replication axis (independent CNV backgrounds converging) remains untested and is
   arguably the most promising untried direction.
5. **No canonical-locus check was completed** against the v3 units, so the
   parallel-adaptation prior (pigmentation, FADS, immune) is untested rather than
   rejected.

## What would overturn this

1. A **parametric or much larger matched null** that reaches genome-scale FDR, showing a
   low-`|X|` locus class at high replication. Our scans cannot see this class even if it
   existed.
2. A **locus class with high between-history heterogeneity (`I^2`) and modest `|X|`**
   surviving matched controls - the better-behaved version of the subset statistic that
   the simulation says is the right target.
3. **Independent (non-nested) histories** with concordant ordered rises: West Eurasian
   turnovers are partly nested by descent, which makes successive rises fewer events
   than they appear. A setting with genuinely independent populations would carry more
   information per observation, and the `1/m!` bound would be more meaningful.
4. **Mutation-level replication** (AMY structural haplotypes, recurrent CNVs) - a
   different kind of independence from drift, and the one direction where the argument
   for repeated directional change is not undermined by shared ancestry.

## Bottom line

**[inference]** The idea that population structure can be used as evolutionary
replication is well-posed, formally sound, and *testable* - and on the best available
data it does not recover a class of loci that Akbari's correction suppresses. The
specific criticism that Akbari "over-corrects" is **unsupported**; the two statistics
agree strongly, and where they could diverge the arithmetic of independent histories
(4-6, not 13) and the overlap of the two blind spots explain why.

The durable products are: a formalized and Lean-verified statement of what the
correction removes, a validated simulation that identifies exactly one regime of
advantage and prices its false-positive cost, an honest three-regime empirical test on
the decisive data, and a corrected reframing (heterogeneity detection versus trend
detection) that any follow-up should start from.
