# A. Scientific summary

> Owner tasks: `recon-akbari`, `recon-lit`, `synth-assessment`. Status: in progress.
> Claims are labelled **[observation]** (measured/derived from data or a proved
> formal result), **[inference]** (a conclusion that follows from stated
> assumptions), and **[speculation]** (a hypothesis not yet tested in this
> project). Numbers that depend on data acquisition or downstream runs are
> written `TBD` with the task that will produce them.

## The question

When several genetically distinct human populations are followed through time,
does the same genetic variant repeatedly change in the same direction along more
than one of those histories? If it does, does that repeated change carry
information that a standard population-structure-corrected selection scan
discards?

This project exists because the current large scan of West Eurasian ancient DNA,
Akbari et al. (2026), corrects for population structure and then asks whether a
variant's allele frequency tracks time. That correction is necessary: migration
and the replacement of one population by another move allele frequencies with no
selection involved. But the same correction can also attenuate a genuine signal
whenever selection and population history are correlated. We build and stress-test
an alternative that uses structure as *replication* rather than only as nuisance,
and then compare it with the structure-corrected scan on the same loci.

## Terms on first use

- **Allele frequency** — the proportion of gene copies at a locus in a sample that
  carry a particular variant; written `p`, in `[0, 1]`.
- **Locus / variant** — a position in the genome, here usually a single-nucleotide
  polymorphism (SNP) or a short insertion/deletion (indel).
- **Population trajectory** ("trajectory") — a set of dated samples that are
  interpreted as one continuous population history, i.e. the population's ancestry
  persists across the sample dates. A trajectory is an **inference**, not a given
  archaeological label. "Inferring trajectories" is the problem of grouping and
  ordering dated individuals into a small number of such histories from genotype
  data alone.
- **Allele-frequency path** — the sequence of allele frequencies of one locus
  along one trajectory's dates.
- **Directional change** — a nonzero, consistently signed trend of allele frequency
  through time (rising or falling), as opposed to a flat or directionless path.
- **Selection** — a systematic, locus-specific fitness difference that drives an
  allele frequency in a consistent direction. Contrast **drift** (random change
  from finite population size), **migration / gene flow** (movement of individuals
  or gene copies between populations), and **admixture/turnover**.
- **Population structure** — the existence of genetically distinguishable groups
  within a sample, produced by ancestry, geography, or drift.
- **Demographic replacement / turnover** — one ancestry stream replacing another
  in a location, so the pooled (mixture) sample's allele frequencies change even
  though neither stream's own frequency need change.
- **Structure correction** — any procedure that removes or downweights
  allele-frequency variation attributable to structure, e.g. ancestry covariate
  adjustment, principal components (PCs), or a genetic relatedness matrix (GRM)
  entered as a random effect.
- **GLMM** — generalized linear mixed model: a regression with fixed effects (here
  a time trend) plus random effects (here ancestry). Akbari et al. fit a binomial
  GLMM with a logit link.
- **Concordance** — the fraction of informative trajectories whose estimated trends
  share the majority direction. A tie counts as one half; exactly-flat slopes are
  excluded.
- **Effective number of independent trajectories, `N_eff`** — a discount of the raw
  trajectory count for shared ancestry/relatedness, so near-duplicate histories do
  not masquerade as independent replication.
- **Replication (in this project's sense)** — the same-direction change observed
  along several trajectories that are **not** close relatives of one another.

## Why a time trend is ambiguous

Naively, a variant whose frequency rises smoothly over time looks selected.
**[observation]** A pooled frequency path can rise monotonically from turnover
alone: if population A (frequency `p_A`) is progressively replaced by population B
(frequency `p_B`), the mixture frequency `(1 − w(t)) p_A + w(t) p_B` moves with the
mixing weight `w(t)` while both source paths are exactly constant. The Lean model
`ReplacementOnly` states this and proves the within-source paths constant
(`lean/TrajectorySelection/Defs.lean`, `Replacement.lean`). So a single pooled
trend is not by itself evidence of selection — the project's guardrails require
more (within-trajectory change, repeated cross-trajectory direction, excess beyond
mixture, independent mutations to the same functional state, or independent
geographic/cultural replicates).

## The two scenarios

**Scenario A — demographic replacement.** Allele frequency changes because the
*mix* of populations changes. The variant may sit at different frequencies in two
source populations that never change their own frequency. A scan that pools the
samples sees a trend; a scan that follows each source separately sees none.
**[observation]** This is a model-level fact (above). Whether a given real locus
is "Scenario A" or "Scenario B" is an **[inference]** that depends on correctly
inferred trajectories and is exactly what the method must decide.

**Scenario B — replicated directional change.** Several genetically distinct
trajectories each move the *same* allele in the *same* direction. Because these are
distinct histories, a single replacement event cannot explain all of them at once.
Repeated movement along independent histories is therefore stronger evidence than
one broad trend. **[inference]** A minimal formal statement (`Replacement.lean`)
shows that, *given* an explicit sign-detection model in which true positive drift
raises the probability of observing a `+` sign above one half, the probability of
perfect two-trajectory concordance strictly exceeds the fair-sign null. That
"given" is essential: the model does not derive the sign bias from physics, and it
covers only two trajectories, independent errors, and sign-only (not magnitude)
evidence.

## Why a structure correction can remove real selection

If a variant's frequency is driven by selection that is *correlated with ancestry*
(for example, selection in an incoming population, or on an allele that tracks a
lifestyle or environmental shift that co-occurs with migration), then the time
trend and the ancestry random effect are partially confounded. Residualizing time
against ancestry — the essence of a GRM-based GLMM — then removes part of the very
signal, and the fixed time-slope variance is inflated.

**[observation]** For a one-mode working covariance `V = σ²I + λuuᵀ` with time `t`
written as `t = αu + r` (`uᵀr = 0`), the finite GLS information for the time slope is
`tᵀV⁻¹t = ‖r‖²/σ² + α²/(σ²+λ)`. The information fraction relative to isotropic
noise is `(1 − cos²θ) + cos²θ · σ²/(σ²+λ)`, where `cos²θ = α²/‖t‖²`
(`lean/TrajectorySelection/Collinearity.lean`). This says the *geometry* of
alignment and the *size of the ancestry mode relative to noise* must **both** be
extreme to erase the slope.

**[inference]** Near-alignment with a dominant ancestry mode therefore predicts
imprecise, under-powered time slopes. **[speculation, to be tested]** That this
materially suppresses *real* Akbari loci is the project's central hypothesis; it is
**not** established by the algebra and must be tested by simulation and empirical
controls (Decision 3, task `sim-sweeps`; Decision 5, task `locus-akbari`). The Lean
result is a precision statement about a linearized working model, not proof of bias
or of a deleted signal in the real fit.

## What would change our minds

We will judge the hypothesis against pre-stated falsifiers rather than confirm it
by construction:

1. **Trajectories not identifiable.** If inference from genome-wide variation
   (excluding the tested locus) does not recover stable, genealogically separated
   trajectories, the design is not viable. (Decision 1, `traj-validate`.)
2. **Uncalibrated null.** If, under neutral demographic simulations and matched
   genomic control loci, the replication statistic does not produce a calibrated
   false-positive rate, we fix the statistic before touching biology. (Decision 2,
   `sim-sweeps`, `sim-ctrl-null`.)
3. **No power separation.** If the method does not recover simulated repeated
   selection better than the Akbari GLMM without a large false-positive increase,
   then it has not solved the problem. (Decision 3, `sim-sweeps`.)
4. **Insensitive controls.** If known loci (AMY, ~20–50 canonical loci) do not
   behave sensibly, we investigate before scaling. (Decision 4, `locus-amy`,
   `locus-canonical`.)
5. **No over-correction class.** If the Akbari reanalysis reveals no meaningful
   low-Akbari/high-replication locus class, we conclude the hypothesised
   over-correction is limited under this test. (Decision 5, `locus-akbari`.)

An honest "unsupported" outcome is a valid result (`docs/MISSION.md`, §Guardrails).

## What the project does not claim, and how it relates to existing work

We do **not** claim that population expansion is selection, that a monotone path is
selected, or that Akbari et al. over-correct without controls. We do **not** treat
archaeological labels as ground truth. **[observation]** A literature survey
(`docs/B-literature-map.md`) found no prior method that jointly infers unlabeled
temporal trajectories, estimates per-trajectory directional change, tests
same-direction change across trajectories with an explicit non-independence
correction, and quantifies what a pooled structure correction suppresses relative
to that replication evidence. **[inference]** The proposed combination is novel as
a design even though each ingredient (trajectory inference, phylogenetic/relatedness
meta-analysis, convergence statistics, structure-aware nulls) has a static or
single-population precedent. **[speculation]** Whether this design in fact finds a
locus class that the Akbari scan misses is the open empirical question the rest of
the document series (C–J) addresses.

## Status of results

All quantitative results are produced by the code inventoried in
`docs/F-reproducible-code.md` and reported, with seeds and provenance, in
`docs/G-validation-report.md`, `docs/H-amy-analysis.md`,
`docs/I-akbari-reanalysis.md` and the decision register `docs/DECISIONS.md`.
Files that do not yet have results carry `TBD` and name the responsible task. No
number in this document is a measurement; the only concrete results here are the
proved algebraic statements in `lean/` and the model-level constructions above.
