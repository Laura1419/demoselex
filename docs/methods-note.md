# Methods note: ARG-based demographic control, and how to apply it downstream

Purpose: state plainly what the method computes, what it has been shown to do on controls,
what it must not be used for, and — the point of the exercise — **how to apply it to the
256 segmental-duplication / CNV loci**.

---

## 1. The question

A structure-corrected association model (GLMM / mixed model with a genetic-relatedness
matrix) removes demographic covariance by projecting it out. That is efficient, and it also
removes signal: an allele that rose because its carrier population expanded is
indistinguishable from one that rose because it was selected, because both appear as
frequency change on the axis the model is correcting.

The alternative posed here is **conditioning rather than projection**: describe the local
genealogy, and ask whether an allele's frequency changed *within* a genealogical background
rather than across backgrounds. Carriage (demography) changes which background people carry;
selection changes allele frequency inside a background.

## 2. The method, in plain terms

1. **Build a dated genealogy.** Take phased haplotypes (here: 3,202 people from the 1000G
   high-coverage release -> 6,404 haplotypes), infer the ancestral recombination graph with
   `tsinfer`, and put times on its nodes with `tsdate`. Genome-wide, per chromosome, in
   10 Mb chunks.
2. **Per variant, measure the age of its carriers.** The variant sits on a branch; the
   samples below that branch are its carriers; the time of the branch's lower node is when
   their common ancestor lived.
3. **Standardise against comparable clades.** Compare that age with clades of the *same
   size* (number of descendant samples), in the *same local recombination regime*
   (breakpoints within +/-250 kb), in the *same compartment* (segmental duplication or not).
   The result is a percentile `p_age`, not a p-value.
4. **Aggregate per window** (1 Mb): count variants in the genome-wide bottom 5% of `p_age`
   and compare with expectation; FDR across windows.
5. **Condition on population.** The same analysis restricted to a labelled population (or,
   better, to ARG-derived local backgrounds, section 4).

## 3. What it shows on controls (EUR, 1,266 haplotypes)

| locus | frequency | clade age | variant percentile | 1 Mb window | reading |
|---|---|---|---|---|---|
| LCT rs4988235 | 0.54 | ~500 gen | 3.8% | 1.48x, q=0.0008 | recent in-situ rise |
| HERC2/OCA2 rs12913832 | 0.18 | 37 gen | 2.0% | q=0.48 | recent rise (young clade) |
| SLC24A5 rs1426654 | 0.57 | ~8,300 gen | 89.7% | q=0.15 | old allele, merely common |
| SLC45A2 rs16891982 | 0.93 | ~2,300 gen | ~90% | 0.99x, q=1.0 | old allele |
| FADS rs174537 | 0.30 | old clade | - | 1.06x, q=0.75 | old allele |
| **AMY cluster** (chr1:103.60-103.72 Mb) | - | youngest 7.6 gen | 5/11 below 5th pct | **1.39x, q=0.0017** | CNV-haplotype structure |

Two things to note. First, the **regime split works**: recent in-situ rises flag, old
alleles that are merely common do not. Second, **population conditioning is essential**: the
identical analysis on the pooled global sample gives LCT 0.9x (p=0.83), because pooling
dilutes its frequency, deepens its carrier clade and fills the comparison set with young
population-specific variants.

## 4. Local backgrounds — the label-free part (the novel piece)

For a locus, each haplotype's **local background** is the most recent clade in the local
tree containing it whose root is younger than a chosen depth tau (walk up from the haplotype
until a node older than tau). No population labels are used.

- The clusters recover real structure on their own (at LCT: an Africa-dominated cluster, an
  East-Asia-dominated cluster, a European + admixed-American cluster).
- **LCT is fixed (frequency 1.00) inside a single 895-haplotype background**, and absent from
  every other large background.
- The **breadth statistic** (how many backgrounds carry >=5% of carriers) separates regimes:
  LCT 1, HERC2 1, SLC45A2 2, FADS 3, **SLC24A5 6**. Breadth is stable across tau; the
  dominant-share numbers are not.

## 5. Demographic baseline

Ne(t) from the ARG itself (generalised skyline, `Ne = k(k-1)dt / 4c`): 126 (0-20 gen), 1,158
(20-100), 1,884 (100-500), 2,779 (500-2,000), 6,293 (2,000-10,000). This is the demographic
context in which any per-locus age should be read, and it is estimated from the same data
rather than assumed.

## 6. Limits — do not claim past these

1. `p_age` is an empirical percentile, not a p-value; windows are not independent (linkage).
2. Absolute ages are substrate-dependent (LCT: 193 gen on a smaller imputed subset, ~500 gen
   here). Percentiles and window enrichments are the robust statements.
3. The window statistic is fragile in low-recombination / pericentromeric regions; the
   per-variant percentile is the reliable instrument.
4. **Allele orientation**: the build passed placeholder alleles to `tsinfer`, so its
   "derived" state is whichever allele is the *global* minority. For a derived allele that is
   the majority in the analysed population, the measured lineage is the complement. This must
   be fixed (carry REF/ALT through the build) before any locus-specific claim about a
   near-fixed allele.
5. Segmental duplications: paralogous variation against a linear reference mimics young
   clades. A CN state must be *genotyped*, not inferred from paralogous SNPs.
6. A related published test (Relate: carriers given the number of lineages at the origin,
   analytic null) was also implemented; it behaves correctly for three of four controls but
   inherits the same orientation issue and needs k0-stratified calibration before genome-wide
   use. It is not required for the downstream plan below.

## 7. Applying it to the 256 segdup / CNV loci

This is the reason the framework exists.

**Inputs.** Per-locus copy-number genotype for every sample (from your genotyping), and a
dated ARG for the same cohort (built from SNPs; the CN states are added as non-inference,
multiallelic sites and overlaid on the topology). If the CNV cohort and the ARG cohort
differ, use the ARG to define local backgrounds for the ARG samples and take background
composition from the shared samples.

**Per locus, compute four things:**

1. **Breadth** — how many ARG-derived local backgrounds carry the CN state (>=5% of carriers).
   One background = a haplotype-confined state; many = a state distributed across backgrounds.
2. **Within-background frequency** — is the CN state fixed inside its dominant background, or
   present at intermediate frequency?
3. **Carrier-clade age percentile** — how young is the carrier clade compared with clades of
   the same size in similar recombination contexts (stratify by compartment = segdup).
4. **Population specificity** — repeat 1-3 within each labelled population; disagreement
   between the label-based and the local-background results is itself informative (it means
   ancestry at that locus differs from the label's assumption).

**Decision rule per locus:**

| pattern | reading |
|---|---|
| young clade + confined to one background + high within-background frequency | candidate recent selection on that CN haplotype |
| old clade + broad across backgrounds | structure / carriage |
| young clade + broad across backgrounds | recent expansion of the CN state across ancestries (interesting; check for artefacts) |
| old clade + confined | old, ancestrally restricted haplotype (drift or ancient selection) |

**Pitfalls specific to this application.** Genotype the CN from reads/graphs, never from
paralogous SNPs; check mapping quality and HWE-like diagnostics in the duplicated region;
report the locus's recombination context (the segdup compartment null); and state the
orientation restriction (limit 4) if any locus's CN state is near-fixed.

**What to present.** A per-locus row: breadth, within-background frequency, clade-age
percentile, population specificity, and the Ne(t) context for dating. Plus the control table
(section 3) as calibration, and the skyline as the demographic frame.

## 8. Provenance

Borrowed: `tsinfer` / `tsdate` (Kelleher et al.; Wohns et al.); the allele-age framing
(Slatkin & Rannala 2000; Griffiths & Tavare 1998); Relate's selection test (carriers given
lineages at origin, analytic null; Speidel et al. 2019); for context, CLUES (Stern, Wilton &
Nielsen 2019) and SDS (Field et al. 2016). Note: **Relate has no omega statistic**; omega /
OmegaPlus is an unrelated haplotype-based method.

Ours, and therefore to be described as an implementation rather than a new method: matching
clades on size x local recombination x compartment; the compartment split; the local
background definition (most recent clade under tau) and the breadth statistic; population
stratification; the window FDR aggregation.
