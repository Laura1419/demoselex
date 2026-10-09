# D. Proposed method

> Owner tasks: `lean-*` and `sim-stat` (equations), `traj-prototype`/`traj-alt`
> (Stage A), `sim-ctrl-null`/`sim-sweeps` (calibration). Status: implementation and
> formalization in progress. Every executable equation below is stated to match the
> code in `sim/trajectory_stat.py`, `sim/empirical_null.py` and `sim/akbari_glmm.py`.
> Formal statements and proofs are in `lean/TrajectorySelection/`. Claims are
> labelled **[observation]** (proved or directly computed), **[inference]** (follows
> from stated assumptions), and **[speculation]** (design hypothesis, not yet
> validated). Unknown quantities are `TBD` with the owning task.

The method has two stages. **Stage A** infers a small set of temporal population
trajectories from genome-wide variation that *excludes the locus being tested*.
**Stage B** estimates allele-frequency change along each trajectory and asks whether
the change is directionally replicated across trajectories, discounting for the
trajectories' shared ancestry. Stage B is implemented now; Stage A is designed and
owned by tasks that are not yet complete, so its equations are labelled as design.

The overall logic is a hierarchy, **Model 0 → Model 3**: Model 0 ignores trajectory
structure (the naive pooled trend); Model 1 uses per-trajectory slopes and raw
concordance; Model 2 adds the relatedness discount and GLS meta-analysis; Model 3 is
a hierarchical shared-mean/heterogeneity model. The executable code currently
covers Models 0–2; Model 3 is specified but not implemented (`TBD`, tasks
`sim-sweeps`/`lean-power`). This ordering matches the module docstring and summary
columns of `sim/trajectory_stat.py`; the bullet list in `README.md` is a
higher-level description of the same ladder.

---

## 1. Definitions and notation

**[observation]** The formal definitions are in `lean/TrajectorySelection/Defs.lean`
(imported by `lean/TrajectorySelection.lean`). For an index set `I` and a
population-state type `P`, a **temporal trajectory** consists of a population map
`population : I → P` and a time map `time : I → ℝ`. An **allele-frequency path**
`p : I → ℝ` satisfies `0 ≤ pᵢ ≤ 1` at every index. An index can be a dated sample or
an aggregate time bin. These definitions impose neither ordered/unique times nor
ancestry independence; those are assumptions to state in later models, not
consequences of the data type.

For a finite observed set `S` with integer sample sizes `nᵢ`, write
`n = Σᵢ∈S nᵢ`, `t̄ = Σ nᵢtᵢ / n` and `p̄ = Σ nᵢpᵢ / n`. The sample-size-weighted
least-squares slope is

```
β̂ = Σᵢ nᵢ (tᵢ − t̄)(pᵢ − p̄) / Σᵢ nᵢ (tᵢ − t̄)².
```

Lean returns `none` if `n = 0` or the weighted time variation is zero; it never
silently treats an undefined slope as zero. `directionSign` maps a slope to `−1`,
`0`, or `1` by exact sign.

Index conventions used throughout: `ℓ` is a locus, `t` a trajectory index, `i` a
time bin (distinct observed time) within a trajectory, `g` a relatedness family, and
`m` the number of families. A **positive slope means increasing frequency with
increasing input time**; because ancient sampling times are usually years before
present (BP), a positive slope means higher frequency in *older* samples. All
slopes are **frequency per input time unit**, not log-odds and not a Wright–Fisher
selection coefficient.

---

## 2. Stage A — inferring temporal trajectories (design)

**[speculation/design]** Stage A takes genome-wide genotypes and sample metadata and
returns a small labeled partition of dated individuals into temporal trajectories,
together with (a) each trajectory's time support and (b) a between-trajectory
relatedness/family assignment for Stage B. The defining requirement is
**non-circularity**: trajectory inference must use genome-wide variation
*excluding the tested locus* (and, preferably, its chromosome or a local window), so
the tested allele cannot influence the histories used to test it.

Two inference routes are planned:

- **Primary (`traj-prototype`, `analysis/trajectory_infer.py`).** PCA on genome-wide
  (test-locus-excluded) genotypes, followed by **time-aware clustering**: projection
  of dated samples into a low-dimensional space, then clustering of sample
  trajectories that are smooth in time, allowing a sample's ancestry to differ
  across dates. Outputs a trajectory label per dated sample and a per-trajectory
  time support.
- **Alternative (`traj-alt`, `analysis/trajectory_alt.py`).** A temporal k-nearest-
  neighbor graph, or a Gaussian mixture / clustering within time windows, used as a
  robustness check against the PCA route. Agreement between routes is a validation
  target (`traj-validate`, Decision 1).

Reported diagnostics (**[inference]** targets, values `TBD`): cluster stability
under resampling, sensitivity to the number of components/clusters, and per-cluster
temporal smoothness. **[observation]** A positive-semidefinite covariance kernel
`TrajectoryCovariance` is available in `Defs.lean` to express non-independence, but
no Stage-A theorem uses or estimates it. The finite sign theorem of §5 presumes
**two correctly separated trajectories**; it is *not* an algorithm or guarantee for
finding them.

**Required hand-off from Stage A to Stage B** (all `TBD` until `traj-prototype` and
`traj-relatedness` complete): the trajectory label per time bin; the time value per
bin; and a *prespecified* family assignment `g(t)` (see §4). Inference must be
performed on background markers, never on the tested locus.

---

## 3. Stage B, Model 0 — the unstratified pooled trend

**Model 0** pools all time bins of all trajectories and computes one count-weighted
frequency slope, exactly as in §1 but over the union of bins:

```
β̂_pooled(ℓ) = Σᵢ nᵢ (tᵢ − t̄)(fᵢ − f̄) / Σᵢ nᵢ (tᵢ − t̄)²,
```

where `fᵢ = yᵢ/nᵢ` is the bin's observed allele frequency. In
`sim/trajectory_stat.py` this is the `summary.pooled_slope` column, computed by
`_weighted_slope(data.loc[data.locus == locus], locus)`.

**Plain language.** This is the naive "does the allele change over time?" answer
that ignores which history each sample belongs to. **Model 0 is for comparison
only.** **[observation]** A trend here can be produced entirely by turnover even when
every within-trajectory slope is zero (§8), so a nonzero `pooled_slope` is not
evidence of selection. The Akbari GLMM is, in spirit, a structure-*corrected*
version of this pooled question; Model 0 is the uncorrected version.

---

## 4. Stage B, Model 1 — per-trajectory slopes and concordance

### 4.1 Per-trajectory weighted slope

Input is a long table with columns `locus, trajectory, time, alt_count, total_count`
(`INPUT_COLUMNS` in `sim/trajectory_stat.py`). Observations that share a locus,
trajectory and time are **pooled as counts** first (`_weighted_slope` groups by
`time` and sums `alt_count`, `total_count`). For one locus `ℓ` and trajectory `t`,
with `yᵢ` alternative-allele count, `nᵢ` total count and `fᵢ = yᵢ/nᵢ` in each of the
`n` distinct time bins `i`:

```
t̄    = Σᵢ nᵢ tᵢ / Σᵢ nᵢ                       (weighted mean time)
f̄    = Σᵢ nᵢ fᵢ / Σᵢ nᵢ                       (weighted mean frequency)
Sxx  = Σᵢ nᵢ (tᵢ − t̄)²
β̂_{ℓt} = Σᵢ (tᵢ − t̄) yᵢ / Sxx                  (== Σᵢ nᵢ(tᵢ−t̄)(fᵢ−f̄)/Sxx).
```

The two forms are equal because `Σᵢ nᵢ(tᵢ − t̄) = 0`; the implementation uses the
second (count) form. Requirements enforced by the code: at least two distinct times,
finite data, and integer counts with `0 ≤ alt ≤ total` and `total > 0`; violations
raise rather than returning a fitted zero. `n` in the `slopes` table is the number
of distinct observed time bins and `span` is `max(t) − min(t)`.

**Roundoff guard.** A floating-point cancellation guard sets only cancellation-scale
slopes to exactly zero:

```
ε_round = 32 · ε_mach · max(max_i fᵢ, 1) / (t_max − t_min),
β̂ := 0   if |β̂| ≤ ε_round.
```

This is a numerical guard, **not** a statistical significance threshold.

### 4.2 Working standard error

For independent binomial counts in distinct bins, with the Jeffreys-style plug-in
`p̃ᵢ = (yᵢ + ½)/(nᵢ + 1)` to avoid a spuriously zero SE at boundary counts,

```
SE(β̂_{ℓt})² = Σᵢ nᵢ (tᵢ − t̄)² p̃ᵢ(1 − p̃ᵢ) / Sxx².
```

This **excludes** drift, overdispersion, date uncertainty and ascertainment, and is
a working (Wald-type) sampling error, not a calibrated significance.

### 4.3 Directional concordance

For locus `ℓ`, count over the prespecified QC-passing trajectories:

```
N₊ = #{t : β̂_{ℓt} > 0},   N₋ = #{t : β̂_{ℓt} < 0},   n_inf = N₊ + N₋,
C_ℓ = max(N₊, N₋) / n_inf          if n_inf > 0,  undefined otherwise.
```

Exactly zero or undefined slopes are not informative and are excluded. `C_ℓ ∈
[1/2, 1]`; a tie gives `1/2` and unanimity gives `1`. This is the function
`_concordance` and the `summary.concordance` column; it is the executable form of
the Lean `directionalConcordance`/`concordanceCounts` (`Statistic.lean`). **[observation]**
The Lean count form is proved equal to `Defs.directionalConcordance`, with range and
fixed-denominator monotonicity theorems. **[inference]** Increasing the majority at
fixed informative `n_inf` cannot lower `C_ℓ`; adding a minority while changing
`n_inf` may. **This descriptive statistic alone is neither a p-value nor evidence of
selection.** In `sim/trajectory_stat.py`, `slopes.direction = sign(β̂)` with exact
zero mapped to `0`.

**Plain language.** For each locus, ask how many distinct histories raised the
allele and how many lowered it. Concordance is the winning side's share. A single
history that rises (or a single pooled rise) gives no replication; agreement of
several independent histories does.

---

## 5. Null calibration of the sign statistic

### 5.1 Conditional fair-sign null

**[inference]** Assume the `n_inf` informative signs are jointly uniform over
`{−, +}^{n_inf}` (independent fair signs). Then `K = N₊ ∼ Binomial(n_inf, 1/2)` and

```
C_ℓ = max(K, n_inf − K)/n_inf = 1/2 + |K − n_inf/2| / n_inf.
```

With `b = binom(n_inf − 1, floor((n_inf − 1)/2))`, the exact null mean and variance
are

```
E₀[C | n_inf]    = 1/2 + b / 2^{n_inf},
Var₀(C | n_inf)  = 1/(4 n_inf) − (b / 2^{n_inf})².
```

**[observation]** Lean proves the binomial-weight normalization, the fold identity
`fairNullConcordance = 1/2 + |k − n/2|/n`, the exact finite-binomial mean formula,
and `Var(C) = E(C²) − E(C)²`; focused examples check `n_inf = 2` (mean `3/4`,
variance `1/16`). The compact central-binomial and second-moment identities are
documented derivations, not separate Lean theorems.

**[inference]** Exchangeability **alone** does not give this null: assigning one
shared fair sign to all trajectories is permutation-exchangeable with fair marginals
yet has `C_ℓ = 1` almost surely. The binomial null requires *independent* fair
signs; sign-dependent exclusion, shared history and correlated errors demand matched
simulations or structure-preserving permutation instead. These finite formulas
condition on an informative `n_inf` and are not calibrated p-values under arbitrary
exchangeable demography.

### 5.2 Replacement vs replication (two-trajectory sign model)

**[observation]** `lean/TrajectorySelection/Replacement.lean` defines two minimal
scenarios. **A (replacement):** source frequencies `p_A, p_B ∈ [0,1]` stay constant
within each trajectory while a time-dependent mixing weight `wᵢ ∈ [0,1]` varies, so
the pooled path is `pᵢ = (1 − wᵢ)p_A + wᵢ p_B`; the theorem proves each within-source
path constant (observed nonzero slopes are noisy sign-only estimates, not slopes of
those exact constant paths). **B (replication):** both bounded paths move from `pᵣ`
to `pᵣ + δ` with a common `δ > 0`.

Under A, the model assigns iid fair signs; Lean proves the distribution of `C` under
A equals the pushforward of the four equally weighted sign pairs, hence `P₀(C = 1) =
1/2` (perfect-concordance event). Under B, independent sign detections with common
`q = P(Sᵣ = + | δ) > 1/2` give

```
P_δ(C = 1) = q² + (1 − q)² = 1/2 + 2(q − 1/2)² > 1/2 = P₀(C = 1).
```

**[inference]** The link from positive drift to `q > 1/2` must come from a separate
estimator/noise model or empirical calibration; the theorem does not derive it. The
event is **two-sided** (`q < 1/2` also raises concordance), so it cannot by itself
identify *positive* change. The proof covers two paths, independent errors and
sign-only evidence; it does **not** cover correlated errors, non-exchangeable
demography, misassigned ancestry or false-positive control. Those are prerequisites
for Decision 2 and the power theorem.

---

## 6. Stage B, Model 2 — relatedness discount and GLS meta-analysis

### 6.1 Family collapse

**Model 2** requires a **prespecified** assignment of each trajectory `t` to an
ancestry/relatedness family `g(t)`, fixed *before* counting replicates. Default
(`family=None` in `trajectory_stat`) is one distinct family per trajectory; a real
analysis must supply `family` from background evidence (`traj-relatedness`), never
from the tested locus. Let `G` be the set of distinct families present and `m = |G|`.
Within each family `g`, `sim/trajectory_stat.py` collapses member slope estimates by
an equal-weight mean and assigns the conservative family error

```
b_g = mean{ β̂_{ℓt} : t ∈ g },
s_g = max( max_{t∈g} SE(β̂_{ℓt}),  SD{ β̂_{ℓt} : t ∈ g } ),
```

where `SD` uses the sample standard deviation (`ddof = 1`) and is `0` for a
singleton family. **[observation]** Consequently an exact within-family copy does
**not** increase precision or `N_eff`. This is a sensitivity convention, not a
learned ancestry-clustering estimator.

### 6.2 Effective number of independent trajectories

Let `ρᵢⱼ` be the prespecified correlation of **family-level slope errors**, with
`ρᵢᵢ = 1` and, in the proved theorems, `0 ≤ ρᵢⱼ ≤ 1`. Put

```
R = Σ_{i∈G} Σ_{j∈G} ρᵢⱼ,      N_eff = m² / R   (m > 0; undefined otherwise).
```

In code, `R` is the sum of the supplied `relatedness` matrix restricted to the
present families, and `summary.n_effective = m² / R`. **[observation]** Lean proves
`m ≤ R ≤ m²`, `1 ≤ N_eff ≤ m`, and that increasing all pairwise correlations at
fixed `G` cannot increase `N_eff` (`Statistic.lean`). **[observation] Duplication
invariant:** if an added trajectory is assigned an already represented family, `G`
is unchanged, so `N_eff` is unchanged — Lean's `effectiveTrajectoryCount_duplicate`
— and in particular cannot increase. Similarity *without* a transitive family
assignment is not enough: the naive raw-label design-effect formula can *increase*
on a duplicate. **[inference]** An incorrect family merger conversely undercounts
genuine independent paths; negative correlations violate the `[1, m]` bounds.
`N_eff` is descriptive and must **not** be inserted as a fractional binomial trial
count.

### 6.3 GLS meta-analysis and heterogeneity

The `relatedness` argument to `trajectory_statistic` is a labeled correlation
matrix (pandas `DataFrame` indexed by family). It is validated: symmetric, unit
diagonal, entries in `[0, 1]`, and **positive definite** (Cholesky); a singular
matrix is rejected rather than jittered into pseudo-independence, and exactly
identical families must be merged. Default `R = I` is an explicit working
independence assumption, not established by the data. With family estimates `b`,
family errors `s`, and `V = diag(s) R diag(s)` (i.e. `Vᵢⱼ = sᵢ ρᵢⱼ sⱼ`):

```
μ̂       = (1ᵀ V⁻¹ b) / (1ᵀ V⁻¹ 1)              (inverse-variance GLS mean)
SE(μ̂)   = 1 / sqrt(1ᵀ V⁻¹ 1)
Q       = (b − μ̂ 1)ᵀ V⁻¹ (b − μ̂ 1),   df = m − 1
I²      = max(0, (Q − df)/Q)   if Q > 0, else undefined
z       = μ̂ / SE(μ̂),    p = 2 · Φ̄(|z|)   (two-sided Wald, working)
q_p     = P(χ²_{df} > Q).
```

These are exactly the `meta_slope`, `meta_se`, `q`, `i_squared`, `z_score`,
`p_value`, `q_p_value` columns. `family_concordance = C` recomputed on the collapsed
family estimates `b` (so it discounts copies). With a single family (`df = 0`)
heterogeneity is undefined (`NaN`); concordance is undefined when no informative
signs remain. The code solves the GLS system via a Cholesky factor of `V`.

**[inference]** `V = diag(s) R diag(s)` treats `R` as the correlation of family
slope *errors*; correlation adjustment must use independently estimated background
relatedness and be checked for sensitivity to the family assignment and to
thresholds. **Neither the Wald p-value nor `Q` is neutral-demographic calibration**;
Decision 2 needs matched simulations or structure-preserving controls. The Akbari
GLMM (`sim/akbari_glmm.py`) is the independent comparator against which these
per-locus summaries are contrasted.

---

## 7. Stage B, Model 3 — hierarchical shared-mean and heterogeneity (specified)

**[speculation/design, `TBD`]** Model 3 would replace the fixed-effect GLS mean of
§6.3 with a hierarchical model in which the per-family true slopes are drawn from a
common distribution with a shared mean `μ` and between-trajectory dispersion `τ²`,
so that heterogeneous trajectories are modelled rather than averaged. The planned
deliverables are the executable sweep and its calibration/power assessment in
`sim-sweeps` (`docs/G-validation-report.md`) and the accompanying separation theorem
in `lean-power` (`lean/TrajectorySelection/Power.lean`). No equations are fixed here
until those implementations exist; stating them now would risk the same
implementation–documentation drift this section is written to avoid. The README
ladder ("Model 3 — hierarchical shared-mean + heterogeneity model") refers to this
item.

---

## 8. Null calibration by matched control loci (`sim/empirical_null.py`)

### 8.1 Replication score

Because Models 1–2 yield a directionally signed, family-discounted change estimate,
the matched-control null ranks loci by a **prespecified** score:

```
score(ℓ) = | μ̂(ℓ) / SE(μ̂)(ℓ) | · family_concordance(ℓ)
         = | z(ℓ) | · C_family(ℓ),
```

exactly `abs(z_score) * family_concordance` in `empirical_null`. Higher means a
stronger, directionally concordant, family-discounted change. **[inference]** It is
a ranking score, **not** a likelihood ratio. All candidate and control loci are
fitted by the same `trajectory_statistic` call.

### 8.2 Matching and empirical P

`empirical_null(counts, metadata, candidates, *, family, relatedness, calipers,
exclusion_bp=1_000_000, min_controls=20, max_controls=None)` requires one metadata row
per locus with `locus, chromosome, position` and six QC features
`allele_frequency, recombination_rate, local_ld, imputation_quality, missingness,
temporal_coverage`, all derived **without using the tested outcome**. Matching
enforces, per candidate:

- same `chromosome`;
- `|position − position_candidate| ≥ exclusion_bp` (default 1 Mb, to reduce local
  LD leakage);
- **exact** `(trajectory, time)` support (the same frozenset of bins);
- every numeric feature within a fixed absolute caliper (defaults:
  `allele_frequency 0.05`, `recombination_rate 0.5`, `local_ld 0.1`,
  `imputation_quality 0.05`, `missingness 0.05`, `temporal_coverage 0.05`);
- then a prespecified scaled-L1 distance `Σ_f |x_f − x_f^candidate| / caliper_f`
  selects the `max_controls` nearest controls, ties broken by locus ID.

With `m` controls and control scores, the two-sided upper-tail empirical P is

```
P(ℓ) = (1 + #{control score ≥ score(ℓ)}) / (m + 1),
```

and `percentile` is the candidate's control midrank
`100 · (#{<} + ½#{=}) / m`. **[inference]** The `+1` includes the tested locus in the
exchangeable rank; ties are conservative; the minimum resolution is `1/(m+1)`.
`P` and `percentile` use different tie conventions and need not sum to 1. The whole
`candidates` set is excluded from every other candidate's control pool. Failures
(missing fields, fewer than two families, undefined concordance, insufficient
matches) **fail closed**; time bins must be harmonized beforehand.

**[speculation]** The exchangeability premise is the hard part: it requires
independence and matching of the *whole* neutral distribution, not just the measured
fields, and small `m` limits the tails. The caller must exclude selected, linked,
batch-confounded or otherwise inappropriate controls upstream, and must not relax
matching post hoc based on outcome.

### 8.3 Demographic simulation

**[speculation, validation in progress]** The matched-control null is complemented
by the demographic simulator (`sim/scripts/simulator.py`, six fixed-seed
`sim/scenarios/*.yaml`) which distinguishes locus-specific selection from pure
ancestry turnover (replacement without selection, selection in one population,
shared same-direction selection, opposing selection, selection aligned with
turnover). Values for false-positive rate and power, and the resulting Decision 2/3
verdicts, are `TBD` in `docs/G-validation-report.md` (`sim-sweeps`).

---

## 9. Relationship to the Akbari GLMM, and the collinearity question

**[observation]** Akbari et al. (2026) fit, per locus, a binomial GLMM
(`sim/akbari_glmm.py` reproduces the structure in small samples):

```
y_i | p_i ∼ Binomial(n_i, p_i),      logit(p_i) = β₀ + β_time t_i + u_i,
u ∼ N(0, τ₁ K + τ₂ I),
```

with `K` a sample genetic relatedness matrix excluding the tested chromosome
(leave-one-chromosome-out, LOCO). The implementation uses a penalized quasi-
likelihood (PQL) working Gaussian approximation: working response
`z = η + (y − n p)/w` with `w = n p(1−p)`, working covariance
`H = diag(1/w) + τ₁ K + τ₂ I`, nonnegative variance components optimized by
restricted working-Gaussian likelihood per locus, then fixed effects and BLUP
updated. Time is centered/scaled internally, so `β_time` and its SE are returned
per unit of input time; `z = β_time / SE`. **[inference]** Because the additive time
and ancestry effects act on **log odds**, the algebra below is a local Gaussian/PQL
working approximation used to study *precision*, not an identity for the nonlinear
likelihood.

**[observation]** Residualize time against the intercept and fixed nuisance
covariates, take a unit GRM eigenvector `u` (`uᵀu = 1`), write residualized time
`t = αu + r` with `α = tᵀu`, `uᵀr = 0`, and assume `‖t‖² > 0`. With
`cos²θ = α²/‖t‖²` (`Collinearity.lean`):

```
‖r‖² = ‖t‖² − α² = ‖t‖² (1 − cos²θ)                    (residual_energy_cosine).
```

For the finite one-mode working covariance `V = σ²I + λuuᵀ` (`σ² > 0`, `λ ≥ 0`,
where `λ` is `τ₁` times that GRM eigenvalue), the GLS time information is

```
tᵀ V⁻¹ t = ‖r‖²/σ² + α²/(σ² + λ),
tᵀ V⁻¹ t / (‖t‖²/σ²) = (1 − cos²θ) + cos²θ · σ²/(σ² + λ)   (working_information_fraction).
```

The inverse of this positive fraction is the **finite** working-model variance
inflation factor (VIF). Its hard-projection/high-mode limit is
`VIF_hard = 1/(1 − cos²θ)` (`projectedVIF_eq`), and the near-alignment working bound
(`near_alignment_working_bound`) gives: if `1 − cos²θ ≤ ε` and
`σ²/(σ²+λ) ≤ δ`, then the information fraction is at most `ε + δ`, hence a working
VIF of at least `1/(ε + δ)` (the inversion step is stated here, not a separate Lean
theorem).

**[inference]** Strong finite-model attenuation needs **both** near-parallel geometry
*and* a large mode variance relative to noise; alignment alone is insufficient.
Even at perfect alignment (`cos²θ = 1`), finite `λ` leaves information
`‖t‖²/(σ² + λ) > 0` and a VIF of `1 + λ/σ²` (`perfect_alignment_working`).
**[inference]** A typical full-rank `K` has `col(K)` equal to the whole observation
space, so membership of `t` in `col(K)` alone cannot imply non-identifiability; the
earlier "hard projection deletes the slope" wording is superseded by this finite-
variance statement (see `docs/C-akbari-reconstruction.md`). The formal Lean
definition of `workingInformation` states the spectral information expression but
does **not** prove the matrix-inverse identity, the PQL likelihood approximation, or
empirical loss of power for Akbari loci. **[speculation]** Whether real Akbari loci
are materially over-corrected requires calibrated simulation and empirical contrasts
(Decision 3, Decision 5); the precision algebra alone cannot establish it.

---

## 10. Assumptions, limitations, and what each number means

- Slopes are frequency-per-time, not log-odds and not Wright–Fisher `s`.
- The binomial SE (§4.2) assumes independent sampling per bin and omits drift,
  overdispersion, date error and ascertainment.
- `C_ℓ`, `N_eff` and `score` are descriptive; their p-values are working-Gaussian
  diagnostics.
- Trajectories are replicates only after non-independence is quantified; a PSD
  covariance kernel can express it but does not estimate it.
- `R` and `family` must come from background evidence, and their sensitivity must be
  reported.
- Stage A is a design, not a proved algorithm; Stage B's sign theory presumes
  correctly separated trajectories.
- Model 3 is specified, not implemented (`TBD`).
- All genomic significance claims require the matched-control null and demographic
  simulation, which are themselves validated only in stipulated simulators so far.

---

## 11. One-paragraph plain-language summary

> "We inferred major population trajectories from genome-wide variation independently
> of each tested locus. We then estimated allele-frequency change through time along
> each trajectory and tested whether the direction and magnitude of change were
> concordant across genetically distinct trajectories, while accounting for their
> non-independence. Significance was calibrated against matched genomic control loci
> and demographic simulations."

## Third-pass prioritisation amendment — 2026-10-09 (append only)

[observation] Operator corrected the objective to ancient-locus PRIORITISATION for downstream modern phenotype work; false positives are expected and acceptable. The prior tests above retain their original interpretation, but do not define the current ranking objective. Frozen scoring rule and validation metrics were written BEFORE rankings in docs/prioritization-prereg.md; implementation is analysis/prioritization_score.py.

[inference] New ranking combines a directional-mixture/common-slope LLR with graded consistency, measured-PC N_eff breadth, an explicit node-path ramp/jump shape component and a nonzero-floor carriage weight. The frozen calibrated null is an empirical common-slope variance law, NOT a drift+jump/path likelihood. Component means and HC3 are measured at the locus; direction/amplitude/drift strengths are not optimized there. Four-of-six agreement receives credit, carriage never causes exclusion. Shape remains a descriptive modifier, not proof of sustained selection.

[observation] See docs/prioritization-results.md for mixed positive-control recovery and weak noisy-step separation; no parameter was retuned to improve those outcomes.
