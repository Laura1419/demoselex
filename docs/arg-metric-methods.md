# Current ARG metric — methods report (2026-10-18)

Purpose: a single unambiguous statement of **what the current analysis actually does**,
what it finds, what is borrowed, what is ours, and what is wrong with it — so the next
steps (local backgrounds, ω, population specificity) are extensions and not redefinitions.

## 1. Data and pipeline

| item | value |
|---|---|
| cohort | 1000G high-coverage phased release, 3,202 samples |
| haplotypes | **6,404** (real phase from the source; no coin-flipped heterozygotes) |
| sites | SNVs only, per-chromosome phased VCF, filtered to **MAF >= 5%** |
| sites retained | 102k (chr22) to 576k (chr2); ~6M genome-wide |
| inference | `tsinfer` per chromosome (path compression), `tsdate` inside-outside in 10 Mb chunks, mutation rate 1.25e-8, Ne 1e4 prior |
| coordinates | **GRCh38** |
| cost | inference ~0.08 s/site; chr22 end-to-end ~2 h; chr1/chr2 ~14 h each |
| code | `analysis/g1k_phased.py` (build/date/stats), `analysis/g1k_pop.py` (subsets), `analysis/g1k_null.py` (null/score/region) |

## 2. The statistic (three layers)

**Layer 1 — per variant.** `age` = the tsdate time of the node that carries the
mutation, i.e. the MRCA of the derived-allele carriers in the local tree.

**Layer 2 — standardisation (empirical, within-tree).**
`p_age = (1 + #{comparison clades with age <= observed}) / (1 + #comparison clades)`,
where comparison clades are drawn from the same chromosome's trees in strata:

- `log2(clade size) +/- 1`  (clade size = number of descendant samples)
- `local recombination bin +/- 1`, where local rho = number of tree breakpoints within
  +/-250 kb, expressed per kb
- compartment: segdup (hg38 genomicSuperDups) vs non-segdup

`p_age` is an **empirical percentile, not a p-value**: it has no sampling distribution.

**Layer 3 — regional aggregation.** Per 1 Mb window: `k5` = number of variants whose
`p_age` is below the **genome-wide empirical 5% threshold**; `expected = 0.05 * n`;
enrichment `= k5/expected`; one-sided Poisson p; Benjamini-Hochberg **q** over windows
within a compartment. (`q` here is FDR. Note: `q_unordered` in the ELAI fork is an
unrelated quantity — the unordered diploid ancestry-pair posterior.)

## 3. Population conditioning — currently a GLOBAL LABEL

The EUR analysis **subsets samples by the 1000G superpopulation label** (CEU/TSI/GBR/FIN/
IBS -> 633 samples, 1,266 haplotypes). This label is global and externally assigned: it is
imposed identically at every locus.

It is **not** learned from local ARG structure. That is a known limitation, not a design
feature — see item A under "next".

## 4. Results (calibration, EUR-restricted)

| locus | variant evidence | window evidence |
|---|---|---|
| LCT rs4988235 | 3.8th percentile in its frequency bin | 1.48x, q = 0.0008 |
| HERC2/OCA2 rs12913832 | 37 gen, 2.0th percentile | 1.17x, q = 0.48 |
| SLC24A5 rs1426654 | 8,270 gen, 89.7th percentile | 1.21x, q = 0.15 |
| SLC45A2 rs16891982 | 2,283 gen, ~90th pct | 0.99x, q = 1.00 |
| FADS rs174537 | old clade | 1.06x, q = 0.75 |
| AMY cluster (chr1:103.60-103.72 Mb) | 5/11 variants below 5th pct (8x excess); youngest 7.6 gen | 1.39x, q = 0.0017 |

Control finding: the same LCT analysis on the **global** sample gives 0.9x (p = 0.83).
Population conditioning is what makes the signal appear; global pooling dilutes the
frequency, deepens the carrier clade (carriers span populations) and fills the
frequency-matched null with young population-specific variants.

Regime split: recent-rise loci (LCT, HERC2) flag; old-allele loci (SLC45A2, FADS) do not.

## 5. Provenance — what is borrowed, what is ours

**Borrowed (all of it, for the statistical ingredients):**
- allele age from a coalescent tree: Slatkin & Rannala 2000; Griffiths & Tavare 1998
- ARG-based selection via branch lengths and allele ages: Relate, Speidel et al. 2019 (omega)
- frequency trajectories / selection timing from dated trees: CLUES, Stern, Wilton & Nielsen 2019
- recent selection from modern genomes without an ARG: SDS, Field et al. 2016
- inference machinery: tsinfer (Kelleher et al.), tsdate / unified genealogy (Wohns et al. 2022)

**Ours (composition only):** matching on clade *size* (the choice that made the null
work), matching local recombination, the segdup compartment split, population
stratification, and the window FDR aggregation. This is an implementation, **not a novel
method**, and it is **weaker** than omega: omega uses all branch lengths in time slices;
we collapse the tree to one node's age.

## 6. Known weaknesses (state these in any presentation)

1. `p_age` is a percentile, not a p-value; windows are not independent (LD), and the 5%
   threshold is empirical.
2. Absolute ages are substrate-dependent (LCT: 193 gen on the smaller imputed EUR set,
   495 gen here) — percentiles and window enrichments are the robust statements.
3. The window statistic is fragile in low-rho / pericentromeric regions and still
   surfaces such windows after centromere masking (chr13:10 Mb, chr21:0 Mb, chr14:10 Mb).
4. Population conditioning uses global labels (section 3).
5. Segdup loci: paralogous variation on a linear reference mimics young clades; AMY is
   read as CNV-haplotype structure, not a point sweep, pending genotyped copy number.

## 7. Next (registered, in order)

**A. Local ARG-derived backgrounds.** At each window, compute pairwise shared branch
length (tskit genetic relatedness) among samples, cluster into local backgrounds, and use
*those* as the conditioning variable — learned per locus from the genealogy, not from a
global label. Test: for loci where local ancestry differs from the label, do the two
analyses disagree, and does the local one agree with known ancestry at that locus?

**B. omega reimplementation on tskit.** Estimate Ne(t) from the ARG itself, define
log-spaced time bands, compute observed vs expected branch length per band, aggregate per
locus. Replaces the fragile window statistic and gives a published, comparable statistic.

**C. Population-specificity test** — with A and B in place: is a sweep confined to a
*local* background, and does the background correspond to a labelled population only where
the label is actually correct?

---

## ITEM A — local ARG-derived backgrounds (implemented, 2026-10-18)

`analysis/g1k_local_bg.py`. Definition: for a locus, each haplotype's **local background**
is the most recent clade in the local tree containing it whose root is younger than tau
(walk up from the sample until a node older than tau; that node is the background). No
population labels are used anywhere.

Validation that the clusters are real ancestry (not noise): at LCT the 12 local clusters
(>=20 haps) are AFR-dominated, EAS-dominated, European+admixed-American, etc. - i.e. the
partition recovers known structure, learned locally.

### LCT is confined to a single local background

tau=200 gen: cluster 49000 holds 895 haplotypes (EUR 550 / AMR 191 / SAS 115) and **every
carrier of rs4988235 (1078) is inside it, at frequency 1.00**; all other large clusters are
0.00. The allele is *fixed within one locally-inferred background*.

### Label-free confinement statistic (breadth), and its tau-dependence

| locus | carriers | backgrounds holding >=5% of carriers (tau=100 / 200 / 500) |
|---|---|---|
| LCT | 1078 | 2 / 1 / 1 |
| HERC2 | 1146 | 1 / 1 / 1 |
| SLC45A2 | 1764 | - / 2 / - |
| FADS | 1924 | - / 3 / - |
| SLC24A5 | 3626 | 5 / 6 / 6 |

Recent in-situ sweeps are confined to one background; the old allele (SLC24A5) spreads
across six. Breadth is stable across tau; the dominant-share numbers are not (at tau above
the carrier-clade age the carriers trivially merge into one cluster), so **breadth is the
statistic to use**.

Caveats: cluster granularity depends on tau and on sample composition; breadth is partly
confounded with allele age (old alleles have had longer to spread) - which is the point,
but means breadth should be compared within frequency/age strata, not globally.

---

## ITEM B — omega attempt FAILED (2026-10-18), diagnosis recorded

`analysis/g1k_omega.py` implements a per-branch statistic:
`omega_b = sum_bands [L_{b,j} / E_j(leafbin_b)] - (#bands touched)` with E taken
empirically as the mean raw branch length among branches overlapping (band, leafbin).

**It is NOT Relate's omega and it is degenerate.** Evidence: LCT rs4988235 gives
omega = -1.0 (ratio sum ~ 0), i.e. no signal, and the implementation does not discriminate
the controls.

Two defects:
1. mechanical: the chunk key used `p//10^7` instead of the absolute chunk start, so
   windows were 10 Mb wide starting at position ~0 (HERC2 raised FileNotFoundError(15,2)).
2. conceptual: the expectation is the mean RAW length among branches overlapping a band,
   which penalises a branch that enters a wide band only briefly (LCT's branch spans
   ~495-500 gen inside a 200-500/500-1000 band). Relate's omega uses the coalescent
   expectation for a branch with a given descendant-leaf count within each epoch, which
   accounts for how much of the epoch such a branch should occupy. The "minus bands
   touched" correction is arbitrary.

**Do not use this statistic.** Correct paths, in order of preference:
(a) run Relate on the same cohort and use its omega directly (published, calibrated);
(b) implement Relate's formula faithfully: estimate Ne(t) (generalised skyline from the
    ARG is straightforward: Ne_j = k_j(k_j-1)dt_j / (4 c_j) with k = mean lineages and
    c = coalescences per band), then use the coalescent expectation per (epoch, leaf count);
(c) if staying empirical, the expectation must be a *distribution* of branch lengths for a
    leaf-count class (e.g. quantiles), not a single mean, and the score must not include a
    bands-touched term.

---

## ITEM B — attempt 2/3 (2026-10-18): mechanically fixed, still not discriminative

`analysis/g1k_omega2.py`. Whole-branch statistic: omega_b = L_b / E(leafcount_b), where
L_b = t[parent]-t[node] (zero-length chains collapsed) and E = median over genome-wide
sampled branches with the same log2 descendant-leaf count; per-SNP score = omega_b - 1.

Fixes applied (each was a real bug in attempt 1):
1. chunk key used the chunk INDEX not its absolute start -> windows were malformed.
2. zero-length branches: tsinfer/tsdate leave nodes differing by ~1e-8 generations, so the
   allele's lineage had to be walked past sub-generation branches (otherwise L=0).
3. **null/scan inconsistency**: the null sampled uncollapsed branch lengths while the scan
   used collapsed ones - a ~10x offset that swamped every signal (a neutral window scored
   the same as LCT).
4. recurrent/oriented mutations: naive per-mutation selection picked a singleton recurrence
   at SLC24A5 (leafbin 1, omega 371) instead of the main mutation.

Result after fixes (EUR, 1,266 haps): LCT omega = 1.66 (branch 2,411 gen vs 1,450
expected). But window medians do NOT discriminate:

| window | median score (omega-1) |
|---|---|
| chr2:135 Mb (LCT) | 2.72 |
| chr2:50 Mb (neutral) | 1.92 |
| chr15:27 Mb (HERC2) | 2.97 |
| chr15:47 Mb (SLC24A5) | 3.83 |

Diagnosis: the expectation is stratified on descendant-leaf count ONLY. Branch length
depends on local Ne and local recombination, which differ between windows, so low-Ne
windows inflate omega - the same failure mode that forced rho + size strata in the
clade-age null.

### The demographic ingredient works: generalised skyline Ne(t)

Ne_j = k_j(k_j-1)dt_j/(4 c_j) from the ARG (chr2, EUR, stride 500):

| band (gen) | mean lineages | coalescences | Ne estimate |
|---|---|---|---|
| 0-20 | 1,359 | 73,285 | 126 |
| 20-100 | 736 | 9,338 | 1,158 |
| 100-500 | 349 | 6,434 | 1,884 |
| 500-2,000 | 198 | 5,244 | 2,779 |
| 2,000-10,000 | 107 | 3,628 | 6,293 |
| 10,000-50,000 | 56 | 209 | 146,346 (too few coalescences; unreliable) |

### Registered next step for B

Use the coalescent expectation instead of the empirical median: for a branch with n
descendant leaves in epoch j, E[l] = integral over the epoch of exp(-(k_j-1)t/(2 Ne_j)) dt,
with k_j and Ne_j taken from the LOCAL (per-window) skyline. That is the principled version
and it uses the skyline above. Do not use the current omega for claims.

---

## ITEM B — RESOLVED differently: Relate's ACTUAL test (2026-10-18), and a factual correction

**Correction of record: Relate (Speidel et al. 2019) has NO omega statistic.** I read the
paper: no mention of omega. Omega/OmegaPlus (Alachiotis et al.) is a haplotype-based
composite-likelihood statistic, unrelated. My repeated references to "Relate's omega" were
wrong and are retracted.

Relate's actual selection test (paper, "Powerful tree-based approaches"): condition on
**k0 = the number of lineages present when the mutation arises**, use **n0 = the number of
present-day carriers** as the statistic, with an analytic null robust to population-size
change (p_R).

`analysis/g1k_relate_test.py` implements this: in a dated tree,
`k0 = 1 + #{internal nodes with time > t_node}`, `n0 = descendant samples`, and the
analytic null from the Beta(1, k0-1) law for the carriers' sample fraction:
**p_R = (1 - n0/n)^(k0-1)**. (An earlier empirical null was floor-limited at ~5e-3, which
is exactly where LCT landed; the analytic null has no floor.)

### Results (EUR, 1,266 haplotypes)

| locus | n0 | k0 | p_R | -log10 |
|---|---|---|---|---|
| HERC2/OCA2 rs12913832 | 802 | 24 | 9.4e-11 | 10.0 |
| SLC45A2 rs16891982 | 1182 | 7 | 8.5e-08 | 7.1 |
| LCT rs4988235 | 678 | 16 | 1.0e-05 | 5.0 |
| FADS rs174537 | 440 | 10 | 2.1e-02 | 1.7 |
| SLC24A5 rs1426654 | invalid | - | - | - |

### Two limitations, both now documented

1. **Orientation.** We passed placeholder alleles ('A','T') to tsinfer, so the inferred
   ancestral state is arbitrary. When the derived allele is the majority in the subset
   (rs1426654 is ~98% in Europeans) the ARG's mutation sits on the COMPLEMENT lineage:
   we measure n0 = 3 instead of ~1,240. Fix: carry REF/ALT through the build (store allele
   orientation per site), which requires a re-download of the VCFs.
2. **Window aggregation is still not calibrated**: the analytic p depends strongly on k0,
   and k0 tracks local Ne, so low-Ne windows are inflated (window fractions of p<1e-4:
   LCT 1.23%, neutral chr2:50Mb 0.67%, HERC2 2.46%, SLC45A2 4.07%, FADS 3.73%). The fix is
   to calibrate within k0 strata (or use each window's own k0-matched null), not to compare
   raw p_R across windows.

So B is now **relate's real statistic, implemented, with correct per-variant behaviour on
three of four controls** (the fourth invalidated by orientation). It is not yet a
calibrated genome-wide scanner because of limitation 2.
