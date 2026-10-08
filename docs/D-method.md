# D. Proposed method

> Owner task: lean-* + sim-stat. Status: TODO.

## Definitions

The formal definitions are in `lean/TrajectorySelection/Defs.lean` (imported by
`lean/TrajectorySelection.lean`). For an index set `I` and population-state type `P`,
a temporal trajectory consists of `population : I → P` and `time : I → ℝ`.
An allele-frequency path `p : I → ℝ` has `0 ≤ pᵢ ≤ 1` at every index. The
index can represent dated samples or aggregate time bins. These definitions
alone impose neither ordered/unique times nor ancestry independence; those
are assumptions to state in later models, not consequences of the data type.

For a finite observed set `S` and integer sample sizes `nᵢ`, let
`n = Σᵢ∈S nᵢ`, `t̄ = Σ nᵢtᵢ / n` and `p̄ = Σ nᵢpᵢ / n`. The sample-size-weighted
least-squares slope is

`β̂ = Σ nᵢ(tᵢ − t̄)(pᵢ − p̄) / Σ nᵢ(tᵢ − t̄)²`.

Lean returns `none` if `n = 0` or the weighted time variation is zero; it
never silently treats an undefined slope as zero. `directionSign` maps a
slope to `−1`, `0`, or `1` by exact sign. For a finite set of trajectory slopes,
let `N₊` and `N₋` count strictly positive and negative slopes; directional
concordance is `max(N₊, N₋)/(N₊ + N₋)`, or `none` if both counts vanish.
Zero slopes are excluded, and a tie gives `1/2`. This descriptive statistic
is **not** a calibrated p-value, evidence of selection, or a claim of
independent histories.

In the replacement-only model, two source frequencies `pA,pB ∈ [0,1]`
are constant while mixing weights `wᵢ ∈ [0,1]` vary. The pooled path is
`pᵢ = (1 − wᵢ)pA + wᵢpB`. Thus a pooled slope can arise from turnover
alone. `TrajectoryCovariance` records a symmetric, positive-semidefinite
finite-support covariance kernel on trajectory labels. Its off-diagonal
entries can represent shared history; no normalization, effective replicate
count, likelihood, or covariance estimator is asserted here. Downstream
formal results must add assumptions for null exchangeability, independence,
noise and identifiability explicitly.

The inherited `Core.lean`, `FrequencyPath.lean` and `Concordance.lean`
scaffolds have older convenience definitions and provisional statistical
placeholders; see `lean/AXIOMS.md`. In particular their `check_replication`
is **not** the directional-concordance definition above. No biological
inference should use those placeholders as calibrated statistics.

## Concordance and effective independent trajectory count (`lean-statistic`)

For locus `ℓ`, count `N₊ = #{t : β̂_{ℓt}>0}` and
`N₋ = #{t : β̂_{ℓt}<0}` among prespecified QC-passing trajectories, and let
`n=N₊+N₋`. Exactly zero or undefined slopes are not informative. Define
`C_ℓ=max(N₊,N₋)/n` for `n>0`, undefined otherwise. It lies in `[1/2,1]`;
a tie gives `1/2` and unanimity gives `1`. Increasing the majority *at fixed
informative n* cannot lower `C_ℓ`; adding a minority while changing `n` may.
The Lean count form is proved equal to `Defs.directionalConcordance`, with
range and fixed-denominator monotonicity theorems in `Statistic.lean`.
This descriptive statistic alone is neither a p-value nor evidence of selection.

Assign each observed trajectory `t` an ancestry/relatedness family `g(t)`
**before** counting replicates. Let `G` be the set of distinct assigned
families, `m=|G|`, and `ρᵢⱼ` their pairwise correlation of family-level
change estimates. With `ρᵢᵢ=1` and `0≤ρᵢⱼ≤1`, put
`R=Σᵢ∈G Σⱼ∈G ρᵢⱼ`, `N_eff=m²/R` when `m>0`, undefined otherwise.
For a covariance interpretation additionally require symmetry and positive
semidefiniteness, and verify that estimated genetic relatedness approximates
the covariance relevant to the tested statistic. The variance of an equally
weighted mean of unit-variance family estimates is `R/m²=1/N_eff`.
Lean proves `m≤R≤m²`, `1≤N_eff≤m`, and that increasing pairwise
correlations at *fixed G* cannot increase `N_eff`.

**Duplication invariant:** if a near-identical added trajectory is assigned
an already represented family, the set `G` is unchanged. Thus
`N_eff(G∪{copy})=N_eff(G)`, in particular it cannot increase. Lean proves
`effectiveTrajectoryCount_duplicate` with precisely this family-membership
assumption. Similarity without a transitive family assignment is *not* enough:
the naive raw-label design-effect formula can **increase** on a duplicate.
For example, nine identical trajectories plus one unrelated trajectory give
`100/82`; duplicating the unrelated trajectory gives `144/85 > 100/82`.
Cluster assignment, threshold sensitivity and PSD should be reported; an
incorrect family merger can conversely undercount genuine independent paths.
Negative correlations also invalidate these `[1,m]` bounds. Do not insert
noninteger `N_eff` as an uncalibrated binomial trial count.

**Conditional fair-sign null:** if the `n>0` informative signs are *jointly
uniform* over `{−,+}ⁿ` (iid fair signs), then `K=N₊∼Bin(n,1/2)` and
`C=max(K,n−K)/n=1/2+|K−n/2|/n`. Let
`b=binom(n−1, floor((n−1)/2))`; the exact null mean and variance are

`E₀[C|n]=1/2+b/2ⁿ`,     `Var₀(C|n)=1/(4n)−(b/2ⁿ)²`.

Indeed, `P(K=k)=binom(n,k)/2ⁿ`, and pairing the two tails gives
`Σₖ binom(n,k)|k−n/2|=n b`: for `k>n/2`, telescope
`(2k−n)binom(n,k)=n[binom(n−1,k−1)−binom(n−1,k)]`.
Also `E[(K−n/2)²]=n/4` and squaring removes the absolute value;
subtracting the squared mean yields the variance. Lean proves binomial-weight
normalization, the fold identity, exact finite-binomial mean formula and
`Var(C)=E(C²)−E(C)²`; the compact central-binomial and second-moment
identities here are **documented derivations, not Lean theorems**. Focused
Lean examples check `n=2`: mean `3/4`, variance `1/16`.

Exchangeability **alone** is insufficient for this binomial null: assigning
one common fair sign to all `n` trajectories is permutation-exchangeable
with fair marginals, yet has `C=1` almost surely. Sign-dependent exclusion,
shared history and correlated errors require matched neutral simulations or
structure-preserving permutation. The finite formulas condition on an
informative `n`, assume independent fair signs, and are not calibrated
p-values under arbitrary exchangeable demography.

## Time–GRM collinearity: which variance is downweighted?

Akbari's binomial GLMM has `Yᵢ | pᵢ ∼ Binomial(nᵢ,pᵢ)` and
`logit(pᵢ) = β₀ + β_time tᵢ + zᵢᵀγ + bᵢ`, with a population-identity
random effect `b ∼ N(0, τ₁K + τ₂I)`. Thus the additive time and ancestry
effects act on **log odds**, not literally on frequency. We use a local
Gaussian/PQL working approximation to study *precision*, not an identity
for the nonlinear binomial likelihood. Residualize time against the
intercept and fixed nuisance covariates first. The formal vectors use an
unweighted Euclidean metric; actual PQL weights and other GRM modes would
change the effective metric and require a separate proof.

Take a unit GRM eigenvector `u` (`uᵀu = 1`), write residualized time
`t = αu + r`, `α = tᵀu`, `uᵀr = 0`, and assume `||t||² > 0`.
`cos²θ = α²/||t||²` is the squared time–mode cosine. The Lean theorem
`Collinearity.residual_energy_cosine` in
`lean/TrajectorySelection/Collinearity.lean` proves

`||r||² = ||t||² − α² = ||t||²(1 − cos²θ)`.

If `1 − cos²θ ≤ ε`, `Collinearity.near_alignment` gives
`||r||² ≤ ε||t||²`. In a *hard-projection/high-mode-variance limit*, the
information about a time slope from this residual is `||r||²/σ²` rather
than `||t||²/σ²`; with nonzero residual the corresponding variance
inflation is `VIF_hard = 1/(1 − cos²θ)` (proved by
`Collinearity.projectedVIF_eq`). Perfect alignment makes this limiting
projected information zero. At that boundary the VIF is undefined/infinite;
Lean's division-by-zero convention is **not** interpreted as a finite VIF.

For the **finite** one-mode working covariance `V = σ²I + λuuᵀ` with
`σ² > 0`, `λ ≥ 0` (λ corresponds to `τ₁` times that GRM eigenvalue),
the GLS time information is

`I_time = tᵀV⁻¹t = ||r||²/σ² + α²/(σ² + λ)`.

Relative to an isotropic-noise baseline, Lean's
`Collinearity.working_information_fraction` proves

`I_time / (||t||²/σ²) = (1 − cos²θ) + cos²θ · σ²/(σ² + λ)`.

The inverse of this positive fraction is the *finite* working-model VIF.
Even when `cos²θ = 1`, finite `λ` leaves information
`||t||²/(σ² + λ) > 0` and a VIF of `1 + λ/σ²`; this is
`Collinearity.perfect_alignment_working`. The formal Lean definition of
`workingInformation` specifies the spectral information expression, but
does **not** prove the matrix-inverse identity, the PQL likelihood
approximation, or empirical loss of power for Akbari loci. In a real GLMM,
estimated variance components, sampling, other covariates and eigenmodes
matter. This precision effect alone cannot establish over-correction;
Decision 3 needs calibrated simulation and empirical contrasts.

`lean/TestCollinearity.lean` checks the `(3,4)` time vector against the
first coordinate: `cos²θ=9/25`, residual energy `16`, hard VIF `25/16`.
For perfectly aligned `(3,0)`, `σ²=1`, `λ=8`, finite information is `1`,
not zero. The file compiles without new `sorry`/axioms.

## Stage A: trajectory inference

The finite sign theorem below presumes **two correctly separated trajectories**; it
is not an algorithm or guarantee for inferring them from genomic data. Verify
trajectory stability and genealogical independence before treating fitted slopes
as replicates (Decision 1). A positive-semidefinite covariance kernel in
`Defs.lean` is available to express non-independence, but this theorem does not
use or estimate one.

## Stage B: replicated directional change

`lean/TrajectorySelection/Replacement.lean` defines two deliberately small
scenarios. **A (replacement):** source frequencies `pA,pB ∈ [0,1]` stay
constant within each trajectory while a time-dependent `w(t) ∈ [0,1]` can
change the pooled frequency `(1−w(t))pA+w(t)pB`. The theorem proves each
within-source path is constant; *observed* nonzero slopes are noisy sign-only
estimates, not slopes of those exact constant paths. **B (replication):** both
bounded paths change from `pᵣ` to `pᵣ+δ`, with a common `δ>0`. This implies
positive **true** changes but does not alone imply any particular distribution
of observed estimates.

## Statistics (Model 0 .. Model 3)

The toy statistic is the existing `directionalConcordance` on two nonzero
estimated slopes with signs `S₁,S₂`. Its possible values are `1` (same sign)
and `1/2` (mixed signs). The event called *high* is **perfect concordance**
`C=1`, not a general threshold or a test for positive direction alone. The
model ignores magnitudes, zero slopes, sign ties from an estimator, missing
observations, trajectory-label uncertainty and slope covariance. Its purpose
is to isolate a minimal mathematical separation; it does not replace a
full fitted Model 0–3 or calibrated genomic significance test.

## Null calibration

Assume under A that noisy estimated signs are **independent fair Bernoulli**:
`P(Sᵣ=+)=P(Sᵣ=−)=1/2`. This exchangeable *sign* null is an additional
calibration assumption, **not a consequence of constant paths**. Lean proves
the distribution of `C` in A equals the pushforward of the four equally
weighted sign pairs for *every decidable event* of `Option ℝ`; hence
`P₀(C=1)=1/2`. In B assume independent sign detections with the same
`q=P(Sᵣ=+ | δ)>1/2` (and `q≤1`). This link between positive drift and
upward sign bias must be established by a separate estimator/noise model or
empirical calibration. Lean proves

`Pδ(C=1)=q²+(1−q)²=1/2+2(q−1/2)² > 1/2=P₀(C=1)`.

The event is two-sided: `q<1/2` also produces excess concordance, so without
the assumed sign-detection direction it cannot identify *positive* change.
Nor does the proof cover correlated errors, non-exchangeable demography,
replacement with systematic within-trajectory ancestry misassignment, more
than two paths, or actual false-positive control. Those are prerequisites for
Decision 2 and the downstream power theorem, not results silently inherited
from this toy calculation.
## One-paragraph plain-language summary
> "We inferred major population trajectories from genome-wide variation independently
> of each tested locus. We then estimated allele-frequency change through time along
> each trajectory and tested whether the direction and magnitude of change were
> concordant across genetically distinct trajectories, while accounting for their
> non-independence. Significance was calibrated against matched genomic control loci
> and demographic simulations."
