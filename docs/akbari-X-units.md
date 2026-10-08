# Akbari et al. 2026 — what the released column `X` actually is, and how to define "low Akbari signal"

Task: resolve the units/definition of the released per-variant column `X`
(and companion `S`, `SE`) in the Akbari et al. selection scan, because our loader
found `mean(X) = +0.412` and `62.3%` of variants with `X > 0`, which is not what a
naive "Z-score" reading predicts.

Labels: **[observation]** = read directly out of a source I inspected;
**[inference]** = my reasoning from observations, not stated by the authors;
**[speculation]** = plausible but unsupported; **[not determinable]** = could not
be established from the sources available to me.

Scope note / tooling limitation: this subagent has **no shell**, so I could not run
`gzip -dc` on `data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz`
(a 674 MB gzip; `read` fails with "Cannot create a string longer than 0x1fffffe8
characters"). Everything below therefore comes from the SI PDF text, the preprint
full text, the Dataverse metadata/README, and the authors' PQLseqPy documentation.
The one question that needs the file header is flagged as such in §(e).

---

## Summary

`X` is **not** a raw Z-score. The SI defines it as `X = Z / √CF`, where `Z` is the
naive `s/SE` from the GLMM and `CF` is a globally calibrated scalar (the SI says
`CF = 2.28`, i.e. `√CF = 1.51`); `CF` is *only* a scale correction, it does **not**
centre or sign-correct the statistic. Two things follow for us: (i) because the
divisor is a single positive constant, any genome-wide positive shift in `Z`
survives intact in `X`, and the authors explicitly document that their null
distribution is non-normal, frequency-dependent and biased (imputation
ascertainment + logit link), and that a large fraction of the genome carries real
signal — so a non-zero genome-wide mean of `X` is expected, not paradoxical;
(ii) `|X| > 5.45` is the authors' stated genome-wide-significance threshold
("π > 99%") and `|X| = 3.16` is FDR = 50%, so "low Akbari signal" should be defined
with these **absolute, two-sided** thresholds (or by rank within MAF-matched
strata), **not** by `X ≈ 0` and **not** by the sign of `X`.

---

## Findings

### (a) Precise definition of `X`

1. **Claim (definition).** `X` is the naive test statistic rescaled by an empirically
   calibrated correction factor: `X = Z / √CF`, `CF = 2.28` (so `√CF = 1.51`).
   **Sources:** Supplementary Information §3, p. 44 of `41586_2026_10358_MOESM1_ESM.pdf`.
   **Support:** direct evidence. **Confidence:** high.
   **Verbatim quote:** *"We implement the FDR approach in practice using a newly
   defined statistic X = Z/√CF, where Z is the nominal test statistic for the
   selection signal without full correction for inflation, and CF is an empirically
   calibrated correction factor. We set CF=2.28 (√CF = 1.51) based on this being the
   point at which we estimate that the posterior probability is 99%."*

2. **Claim (what `Z` is).** `Z` is "the number of standard errors this quantity is
   from zero", i.e. `Z = s / SE(s)`, from the per-variant GLMM fit
   `logit(p_ij) = α_j + s_j t_i + MVN(0, σ_j²K)` with `H0: s_j = 0`; `s_j` is
   re-estimated independently at every variant (`σ_j²`, `α_j` free too).
   **Sources:** Main text "Testing for selection while correcting for population
   structure" (preprint `PMC11419161`; same Methods in the Nature version).
   **Support:** direct evidence. **Confidence:** high.
   Note the naming collision: in the PQLseqPy API, `X` is the *covariate matrix*, and
   the package reports `z_beta`, `p_beta` — the paper's `X` is a post-hoc transformed
   `z_beta`. Do not confuse the two when reading code. **[inference]**

3. **Claim (`Z` is not a valid p-value by itself).** The authors state the naive
   statistic is miscalibrated: *"In practice, the statistic needs recalibration as it
   is inflated due to unmodeled features of the data, so we empirically assess
   significance from enrichment of signals in independent GWAS."*
   **Sources:** Methods, preprint `PMC11419161`. **Support:** direct evidence.
   **Confidence:** high.

4. **Claim (what "adjusted for inflation" means).** It means **global rescaling by a
   single empirically calibrated constant**, and explicitly **not** genomic control /
   median-χ² rescaling and **not** the LDSC intercept. The SI evaluates those and
   rejects them: λGC-corrected threshold is "too conservative" (49 independent loci);
   the LDSC "pseudo-intercept" gives CF = 2.33 (statistically similar to 2.28) but
   "because of the non-linearities, we refrain from using this correction to control
   formally for FWER"; the adopted CF comes from the point where the GWAS-enrichment
   curve plateaus (~99% posterior).
   **Sources:** SI §3 "Controlling for family wise error rate (FWER)" and
   "Controlling for false discovery rate (FDR)…", pp. 39–44. **Support:** direct
   evidence. **Confidence:** high.

5. **Contradiction inside the paper (recorded, not resolved):** the main text says
   the threshold occurs at `Z = 9.10` and that the rescaling divisor is **1.67**
   ("we rescale the naïve score by this quantity to obtain an X-statistic (Z/1.67)"),
   whereas the SI says `√CF = 1.51` (`CF = 2.28`). `9.10/5.45 = 1.67`; `5.45 × 1.51 =
   8.23`. The SI's own Figure S3.2 caption is consistent with **1.51** ("The black
   dashed line shows the P-threshold (1.95e-16) with a posterior probability of 99%";
   two-sided p = 1.95e-16 ⇔ |Z| ≈ 8.23), and so is the referees' summary in the
   peer-review file ("the empirical threshold (|X| > 5.45, corresponding to P-value
   < 1.95 × 10-16)"). Either way, **the authors' stated operating threshold is
   `|X| > 5.45`**, which is what our classification needs. **Sources:** main text
   (preprint `PMC11419161`); SI §3 p. 44 and Fig. S3.2; Springer Nature peer-review
   file `41586_2026_10358_MOESM4_ESM.pdf` (quoted in search results).
   **Support:** direct evidence for both statements, **[inference]** that the released
   file follows the SI (1.51). **Confidence:** high that the two texts disagree;
   medium on which one was used to generate the released numbers.

### (b) Is `X` expected to be centred at zero? (No — and why)

6. **Claim.** The authors never state that `X` is centred at zero, and they describe
   three properties that make a mean-zero expectation wrong:
   (i) *"In this study, however, the null distribution is not normal … making FWER
   inappropriate for defining significance."*
   (ii) *"a substantial fraction of the genome appears to be affected by real
   directional selection (Extended Data Figure 2), violating the assumption that most
   variants are null."* — quantified as **45% of the genome** in at least weak LD
   (r² > 0.05) with <5%-FDR signals, and **2.35 ± 0.13%** of the logit-scale allele
   frequency variance attributable to directional selection.
   (iii) Frequency-dependent **bias** (not merely variance inflation) in the
   Z-scores: *"GLMM summary statistics exhibit frequency-dependent bias and inflation
   in Z-scores. A major contributor is SNP ascertainment due to imputation using
   modern reference panels, which preferentially retain variants with higher minor
   allele frequency in the reference panel at time zero (MAF0) … SNPs that drifted
   upward in frequency are more likely to pass quality control, creating an apparent
   bias and inflation in Z-scores … this ascertainment mechanism alone can produce a
   genomic inflation factor even in the absence of population structure or true
   selection."* Figure S2.13 is captioned *"The left panels show mean Z-scores and
   the right panels show genomic control inflation factors (λGC) as a function of
   minor allele frequency"*, i.e. the authors themselves plot a frequency-dependent
   **mean** Z, not just a median.
   **Sources:** SI §2 "Non-normality of the null distribution in GLMM test statistics"
   (pp. 34–35) and Fig. S2.13; SI §3 "Controlling for FDR…" (p. 41); main text.
   **Support:** direct evidence for the quotes; the following is **[inference]**:
   a positive genome-wide mean shift is exactly the direction implied by
   "ascertainment of variants that drifted upward", and a global division by √CF
   cannot remove it. **Confidence:** high for the mechanism and for "not centred at
   zero in general"; **medium** that it quantitatively explains our `mean(X) = +0.412`
   / 62.3% positive (see §Contradictions and §Missing evidence).

7. **Our own observation, not verified by me.** `mean(X) = +0.412`, `62.3%` of
   `X > 0`, over 9,739,624 variants — supplied by the parent's loader.
   **[observation — parent loader, unverified in this subagent]**. The release does
   contain exactly 9,739,624 variants, matching SI §1: *"A total of 9,739,624
   variants, including 8,074,573 SNPs and 1,665,051 indels, passed all three steps of
   QC"* (**direct evidence**, SI §1 p. 7), so the row set is right.
   Interpretation **[inference]**: the offset is a mixture of (i) the ascertainment /
   logit-link bias the authors document, (ii) real pervasive directional selection
   (their own claim that selection affects an overwhelming majority of the genome),
   and (iii) possibly a column or loader mismatch on our side, which the file header
   would settle (see (e)). Note a `+0.412` mean on `X = Z/1.51` corresponds to a mean
   `Z ≈ 0.62`, which is large, so option (iii) deserves a 5-minute check.

### (c) Exact significance thresholds; what `|X| > 5.45` is

8. **Claim (thresholds).**
   - `|X| > 5.45` ⇔ **π > 99%** posterior probability of selection ⇔ the study's
     genome-wide-significance threshold. *"We identified 347 independent loci (279
     excluding the HLA region) with |X|>5.45, corresponding to a π>99% probability of
     selection."* The larger, final list uses `|X|>5.45` too: 479 independent loci
     (410 non-HLA) with π>99% (SI §5).
   - `|X| = 3.16` ⇔ **FDR = 50%**: *"Using a threshold of |X|=3.16, which
     corresponding to FDR=50%, we identify 10361 non-HLA loci."*
   - Binary annotation used for the FDR machinery: SNPs with FDR < 1% vs ≥ 1%
     (Fig. 1c legend).
   **Sources:** main text (preprint `PMC11419161`); SI §5 (p. 47).
   **Support:** direct evidence. **Confidence:** high.

9. **Claim (what `|X|>5.45` means numerically).** The authors *intend* `X` to be
   interpretable against a standard normal: *"With our parameterization of the
   X-statistic, the threshold |X|=5.45 which if it was a normally distributed
   variable would imply a classically genome-wide significance threshold of 5x10-8,
   is the threshold of genome-wide significance for our study as well."*
   **But note** that with the SI's own `√CF = 1.51`, `|X| = 5.45` corresponds to
   `|Z| = 8.23`, i.e. a **two-sided p ≈ 1.95 × 10⁻¹⁶**, not 5 × 10⁻⁸ (this is exactly
   the value in Fig. S3.2 and in the referee summary). The citation-friendly
   statement is therefore: *"|X|>5.45 is the paper's genome-wide-significance
   threshold, with π>99%; the paper equates it with a normal-theory p-value at the
   GWAS convention."* **Sources:** SI §3 p. 44; Fig. S3.2; peer-review file.
   **Support:** direct evidence (both numbers) + **[inference]** on the arithmetic.
   **Confidence:** high.

### (d) Units and sign convention of `S` (the selection coefficient) and `SE`

10. **Claim (sign convention).** `s` is the coefficient on time for the **alternative
    (ALT) allele** — not the derived allele, not the ancestral allele, not "the
    reference allele". The GLMM models the ALT allele as the "success":
    *"Y[i, 0] = number of **alternative alleles**; Y[i, 1] = number of **reference
    alleles**"* and `logit(p_i) = X[i]β + u_i` where `p_i` is the ALT frequency
    (PQLseqPy README, mokar01/PQLseqPy, the frozen version of which is deposited as
    `PQLseqPy-v0.1.2.zip` in the Dataverse release). The release's variant IDs are
    `CHROM_POS_REF_ALT`, so the oriented allele for every row is recoverable from the
    ID. **Sources:** PQLseqPy README (mokar2001/PQLseqPy); Dataverse README.txt
    ("Variant IDs are defined as CHROM_POS_REF_ALT"). **Support:** direct evidence.
    **Confidence:** high. Caveat **[inference]**: the *paper figures* label their
    trajectories with "derived allele frequency" and `Anc=<base>` (e.g. SI §5 panels),
    which is a **different orientation** from the regression's ALT orientation; the
    methods note the two panels' derived-allele labelling should not be read as the
    definition of `sign(s)`.
11. **Claim (time unit / what a positive `s` means).** `t_i` is *"the negative
    sampling date in the past, in units of twice the generation interval (assuming 29
    years per generation)"*; `s_j` is *"the per-generation selection strength at the
    allele, assumed constant over time and space"*, so `s > 0` = ALT allele log-odds
    **increasing toward the present** (present-day individuals have `t = 0`, ancient
    individuals `t < 0`). **Sources:** Methods, preprint `PMC11419161`.
    **Support:** direct evidence for the quotes. **Confidence:** medium — see
    Contradictions #2: "per-generation" and "units of twice the generation interval"
    are not mutually consistent at face value, and I could not resolve it from the
    text I could reach. Practical check: the reported median `|s| = 0.8%` (range
    0.4–4.2%) with `SE ≈ 0.1%` for common variants is the magnitude expected for a
    **per-generation** coefficient, so the printed numbers are most plausibly
    per-generation **[inference]** (the Nature version reports `SE` "typically 0.17%
    for common variants" in the same sentence).
12. **Claim (`SE`).** `SE` is the standard error of the estimated `s` (the same
    quantity whose ratio defines `Z`), with *"a standard error typically around 0.1%
    for common variants"* (preprint) / *"typically 0.17% for common variants (Extended
    Data Fig. 1c,d)"* (published Nature text). Also reported per-locus in SI §5 as
    `s=0.47%±0.08%`, `s=-0.55%±0.09%`, etc. **Sources:** Methods; published Nature
    article text; SI §5 panels (verbatim example: `s=-0.58%±0.08%, X=-7.07,
    P_X=1.6e-12`). **Support:** direct evidence. **Confidence:** high. I found no
    statement that `s` or `SE` is "scaled by mean frequency"; the logit link makes
    `s` frequency-independent by construction ("The logit function … transforms
    allele frequency so its expected change per generation is proportional to the
    selection coefficient s_j (regardless of p_ij)"), which is a **direct-evidence
    answer to "scaled by mean frequency?" = no** (they are per-unit-time logit
    slopes; not normalised by frequency). [inference on the "so no" step]
13. **Useful calibration anchors for `SE`/units in practice** (verbatim from SI §5):
    `s=0.47%±0.08%, X=5.73, P_X=9.9e-09`; `s=1.25%±0.22%, X=5.78, P_X=7.3e-09`;
    `s=-0.84%±0.09%, X=-9.68, P_X=3.8e-22`. These confirm `X ≈ (s/SE)/1.51`:
    check `0.47/0.08 = 5.9`; `5.9/1.51 = 3.9` (mismatch — see note) and
    `1.25/0.22 = 5.7`; `5.7/1.51 = 3.8` (mismatch). **[inference]** These printed
    `s±SE` pairs give `s/SE ≈ X × 1.0`, i.e. the printed `X` is numerically close to
    the raw `s/SE` for the top loci, not to `s/SE/1.51`. Two non-exclusive readings:
    (i) the printed `±SE` in the figure panels is rounded to 2 decimals and the true
    SE is ~1.5× larger (e.g. `0.47 ± 0.121` → `Z = 3.9` → `X = 2.6`, still not 5.73);
    (ii) the released/plotted `X` may be scaled differently from the SI's stated
    `√CF = 1.51` (consistent with the main text's 1.67 being a *different* quantity
    again). Because the rounding can absorb at most a factor ~1.3, this arithmetic
    suggests the `P_X` values printed alongside (`9.9e-09`) are normal-theory p-values
    of `|X|`-like quantities (`p(5.73) ≈ 1e-08`), i.e. **`X` is internally consistent
    as a normal-theory Z, which itself argues that the applied divisor is near 1, not
    near 1.5.** This is a genuine unresolved inconsistency; flagged in §Missing
    evidence. **Confidence:** medium (it is arithmetic on rounded printed values).
    Consequence for us: **do not assume the released `X` equals `s/SE/1.51`**; use
    `X` only through the authors' own threshold semantics (`|X|>5.45` = π>99%).

### (e) Companion posterior column, and how it relates to `X`

14. **Claim.** The posterior probability of selection `π` is a **deterministic
    function of `X` (its magnitude)** obtained from the FDR/enrichment calibration:
    *"To translate X to a posterior probability of selection π, we used a False
    Discovery Rate (FDR) approach"*, with *"at X-statistic magnitudes greater than our
    threshold for genome-wide significance of 5.45, π>99%"*; Extended Data Figure 2b
    is explicitly *"(b) Mapping X-score to posterior probability (π), False Discovery
    Rate (FDR), number of independent loci excluding the HLA region (N), and the
    percentage of the genome in LD (r² > 0.05)"*; and Figure 3 panel legends print
    `X=…, P_X=…, π>0.99` together. So `π` is a monotone function of `|X|` — it carries
    no information beyond `|X|` apart from the calibration mapping.
    **Sources:** main text; SI §3 pp. 41–44; Extended Data Fig. 2b caption.
    **Support:** direct evidence. **Confidence:** high (as a *mapping*).
15. **Whether `π` is a column in the released file: [not determinable from the
    sources I could inspect].** Neither the SI's Supplementary/Online Table 4
    description (*"Summary statistics for selection in 9.7 million variants. The data
    are provided as a tab-delimited text file, compressed using gzip."*) nor the
    Dataverse `README.txt` lists the columns; the README only says *"Coordinates are
    in GRCh37; Variant IDs are defined as CHROM_POS_REF_ALT"*. Because they are
    redundant with `|X|`, their absence costs us nothing: we can reconstruct
    `π ≈ 0.99` at `|X| = 5.45` and `FDR = 0.5` at `|X| = 3.16`. **What would answer
    it:** `gzip -dc data/raw/akbari/Selection_Summary_Statistics_01OCT2025.tsv.gz |
    head -1` (the header), which I could not run (no shell in this subagent).

---

## Contradictions

1. **Scaling factor 1.51 vs 1.67.** SI §3, p. 44: `X = Z/√CF`, `CF = 2.28`
   (`√CF = 1.51`). Main text: `Z = 9.10` is `1.67`-times `5.45` and *"we rescale the
   naïve score by this quantity to obtain an X-statistic (Z/1.67)"*. Corroborating
   the SI: Fig. S3.2 caption gives the 99%-posterior threshold as p = 1.95e-16
   (⇔ |Z| ≈ 8.23 = 5.45 × 1.51), and the referee report in the peer-review PDF
   describes *"the empirical threshold (|X| > 5.45, corresponding to P-value <
   1.95 × 10-16)"*. Recorded as unresolved; **the operating threshold `|X| > 5.45`
   is the same in both readings**, so our classification is not destabilised.
2. **"per-generation" `s` vs `t` "in units of twice the generation interval".**
   Both in the same Methods paragraph (preprint `PMC11419161`). Unresolved.
3. **Median χ² inflation 5.26-fold (main text) vs 4.75 (SI, "excluding the HLA
   region").** Possibly different variant sets (with/without HLA), but the main text
   does not say "excluding HLA" while the SI does. Recorded.
4. **Printed `X` values vs the SI's stated divisor** (finding 13): for the SI §5
   top-locus panels, `s/SE` ≈ `X` (ratio ~1.0) rather than `X × 1.51`, and the
   printed `P_X` matches a normal-theory two-sided p of `|X|` (`p(5.73) = 9.9e-09`,
   `p(6.26) = 3.8e-10`, `p(9.68) = 3.8e-22` — all reproduce). So whatever constant
   was actually applied to produce `X`, the printed `X` behaves as a standard-normal
   Z. Recorded.
5. **Our loader's `mean(X) = +0.412`, 62.3% positive** vs the naive "Z-score is
   mean-zero under the null" reading of the file caption. The paper documents why a
   mean-zero null fails here (findings 6), but **never states the genome-wide mean or
   sign balance of `Z` or `X`**, and never states that `s > 0` is more common than
   `s < 0`. So the offset is *expected in kind*, not *quantified by the authors*.

---

## Correct implication for OUR analysis

**Q: should "low Akbari signal" mean `X` near zero, or `X` below the genome-wide mean?**

**Neither, as stated — the well-posed definition uses the magnitude `|X|` with the
authors' absolute thresholds, or a rank/quantile within matched strata.** Concretely:

1. **Do not use `X ≈ 0` as the "low" pole.** `X` is `Z/√CF` with a single global
   positive divisor. It is a rescaled test statistic, not a centred effect size, and
   (i) the paper's own null is non-normal and frequency-dependently biased, and
   (ii) our loader sees a large positive offset (`mean = +0.412`, 62.3% positive).
   "Signal" in this statistic is **distance from zero in either direction**, i.e.
   `|X|`, because selection in either direction on the ALT allele gives a large
   `|X|` (all the authors' thresholds and `π` maps are written with `|X|`).
2. **Do not use the sign of `X`.** `sign(X) = sign(s)` = direction of ALT-allele
   frequency change over time. With 62.3% of variants positive, `X < 0` is a
   *selected minority direction*, not "low signal". If we want direction, we must
   also convert alleles: the released IDs are `CHROM_POS_REF_ALT`, so `s > 0` means
   "ALT allele trending up"; flipping to derived/ancestral orientation requires an
   ancestral-state call that this file does not carry.
3. **Do use the authors' own calibration for the "high" pole:**
   `|X| > 5.45` ⇒ `π > 99%` (their genome-wide significant set); `3.16 < |X| < 5.45`
   is "FDR between 50% and 1%" (the `|X| = 3.16` ⇔ FDR = 50% anchor). The "low" pole
   should be **`|X|` far below 3.16**, e.g. the bottom quintile/decile of `|X|`,
   or an explicit band such as `|X| < 1` — choose it **within MAF and
   functional-annotation strata**, because the authors show the statistic's bias and
   variance both depend strongly on MAF (Fig. S2.13: mean Z and λGC plotted against
   MAF; imputation ascertainment depends on MAF0). Unmatched absolute cutoffs will
   confound "low Akbari signal" with "rare variant".
4. **Making the low-Akbari / high-replication quadrant well-posed.** Define the
   Akbari axis as a *strength* axis and make both arms explicitly two-sided and
   matched:
   - `A_high` := `|X| > 5.45` (π > 99%). Optionally an intermediate band
     `3.16–5.45`.
   - `A_low` := `|X|` below an explicitly chosen low quantile (e.g. the lowest 20%)
     **computed within MAF bins** (and, if we care, within the same functional/LD
     context, or with the HLA region excluded/handled separately — the authors
     themselves exclude HLA for locus counting).
   - Replication axis unchanged, but cross-tabulate with MAF and with `SE` (or with
     MAF-matched expected `|X|`), because `|X|` is driven by `s/SE` and `SE`
     depends on frequency.
   - Report sensitivity to: (i) absolute-zero vs median-centred cutoff
     (`|X - median(X)|`), (ii) `|X|` vs signed `X`, (iii) MAF-matched vs unmatched
     thresholds. If the quadrant claim survives all three, it is robust; if it only
     survives with "X near 0" as the low pole, it is an artefact of the calibration
     offset.
5. **A defensible sentence for the paper:** *"Akbari et al. define genome-wide
   significance as |X| > 5.45 (posterior probability of selection π > 99%) and
   |X| = 3.16 as FDR = 50%; X = Z/√CF is a globally rescaled test statistic, not a
   mean-zero Z-score, and the authors document a frequency-dependent, non-normal null
   (and ~45% of the genome in LD with signals), so X is not centred at zero. We
   therefore score 'Akbari selection signal' by |X|, using |X| > 5.45 for the
   high-signal arm and the bottom MAF-matched quintile of |X| for the low-signal arm,
   with the HLA region handled separately."*

---

## Missing evidence

- **The released file's own header/columns.** I could not decompress the .gz (no
  shell). This blocks a definitive statement on (i) whether a `π`/posterior column
  exists, (ii) the exact column names and their order (our loader's `X` column
  identity), (iii) whether `X` in the file is `Z/1.51`, `Z/1.67`, or something else.
  Answer with two commands: `gzip -dc … | head -1` and
  `gzip -dc … | head -100000 | awk -F'\t' '{…}'` for a mean/median/quantile check on
  the actual `X` column, plus a `median(X)`/`mean(X)` comparison against the
  `s/SE` columns if both are present.
- **The MRD of the offset.** No source states the genome-wide mean, median or sign
  balance of `Z`/`X`, nor the MAF-stratified means that Figure S2.13 plots. If we can
  read those numbers off Fig. S2.13 (or re-derive them from the released file
  stratified by MAF), we can decide whether our `+0.412` is the documented
  ascertainment bias or a loader artefact.
- **Exact unit of `t` (and hence of `s`).** "per-generation" vs "twice the
  generation interval" is unresolved from the preprint/SI text; the published Nature
  Methods (PDF at reich.hms.harvard.edu exceeds the 20 MB fetch limit, and
  nature.com's HTML fetch returned only the reference list) would settle it. Also
  unresolved: the printed `s/SE` vs printed `X` arithmetic of finding 13.
- **Whether the `1.51` vs `1.67` discrepancy was corrected in the final Nature
  Methods text.** Only the (inaccessible) published Methods can settle this.

---

## Sources

**Kept (all inspected directly):**
- *Akbari et al. 2026, Supplementary Information* —
  `https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10358-1/MediaObjects/41586_2026_10358_MOESM1_ESM.pdf`
  (also local `data/raw/akbari/41586_2026_10358_MOESM1_ESM.pdf`) — the authoritative
  source for `X = Z/√CF`, `CF = 2.28`, FDR/FWER discussion, the non-normality and
  ascertainment-bias subsection, the threshold statements, and Tables S5.1 / §5
  panels with `s±SE, X, P_X`.
- *Preprint full text, `PMC11419161`* — `https://pmc.ncbi.nlm.nih.gov/articles/PMC11419161/`
  — Methods (GLMM equation, units of `t`, `H0: s = 0`), main-text threshold and
  `Z = 9.10` / `1.67` statement, Online Table 4 description.
- *Nature article landing page* — `https://www.nature.com/articles/s41586-026-10358-1`
  — used only to confirm the published text variant ("SE typically 0.17%"); the
  fetched HTML contained the reference list, not Methods.
- *Springer Nature peer-review file* — `…/41586_2026_10358_MOESM4_ESM.pdf` (quoted via
  search results; direct fetch timed out / "Invalid PDF structure") — independent
  confirmation that `|X| > 5.45` ⇔ p < 1.95e-16.
- *Harvard Dataverse release* — dataset metadata + `README.txt`
  (`https://dataverse.harvard.edu/api/datasets/:persistentId?persistentId=doi:10.7910/DVN/7RVV9N`,
  file id 13596811) — confirms the exact file name/size/categories of
  `Selection_Summary_Statistics_01OCT2025.tsv.gz`, that it is *not* tabular-ingested,
  varID = `CHROM_POS_REF_ALT`, and that the README contains **no** column dictionary.
- *mokar2001/PQLseqPy README* — `https://github.com/mokar2001/PQLseqPy` — the
  ALT-vs-REF orientation ("Y[i,0] = number of alternative alleles") and output
  columns (`beta, se_beta, z_beta, p_beta`), i.e. the sign convention of `s`.

**Rejected / deprioritised:**
- `reikh.hms.harvard.edu` published Nature PDF
  (`…/inline-files/2026_Akbari_Nature_selection_0.pdf`) — exceeds the fetcher's 20 MB
  PDF limit; only usable via search snippets (used for one quote).
- The .xlsx supplementary tables and the 674 MB .gz — no shell available to inspect
  them in this subagent.
- Assorted generic hits (Bartleby z=1.67 homework pages, GWAS-summary-statistics
  column-vocabulary docs, unrelated "adjusted Z-score" accounting papers) — no
  bearing on this paper.

## Next steps (most useful only)

1. **Dump the released header and basic column stats** (5 minutes, parent has a
   shell): confirm the `X` column identity, check for a `π`/posterior column, and
   compute `median(X)`, `mean(X)` and `corr(X, s/SE)` per MAF bin. This closes
   finding 15 and the last open item in Contradictions #4.
2. If a `s` and `SE` column exist alongside `X`, test `X ≈ s/SE/k` for `k ∈ {1, 1.51,
   1.67}` on the top-|X| variants to determine which scaling the release used.
3. Only if the quadrant claim is sensitive to the choice: read Figure S2.13's
   MAF-stratified mean-Z values (or re-derive them) to build a MAF-matched null for
   `|X|`.
