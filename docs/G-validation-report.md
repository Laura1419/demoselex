# G. Validation report

> Owner tasks: sim-sweeps, locus-canonical, traj-validate. Status: IN PROGRESS.

## Simulation Framework (sim-scaffold)

A demographic simulator has been implemented to support validation of the trajectory replication method. The simulator produces dated individuals with genome-wide genotypes along population trajectories with configurable migration, replacement, drift, and selection.

### Features
- **Configurable scenarios**: 6 core scenarios for method validation
- **Provenance tracking**: Records seed, versions, and execution context
- **Output formats**: Dated genotype matrices + truth tables
- **Multiple populations**: Configurable population trajectories

### Core Scenarios
1. **Neutral drift**: Baseline - no selection, only drift and migration
2. **Replacement only**: Population turnover without selection
3. **One-population selection**: Selection in single population only
4. **Shared selection**: Selection shared across multiple populations
5. **Opposing selection**: Selection in opposite directions
6. **Selection aligned with turnover**: Selection in arriving populations

### Technical Implementation
- **Language**: Python 3.12 with NumPy/Pandas
- **Configuration**: YAML-based scenario definitions
- **Genotype representation**: Diploid individuals with binary alleles
- **Selection model**: Additive fitness with configurable strength
- **Migration model**: Fractional individual replacement

### Usage
```bash
python3 sim/scripts/run_simulations.py sim/configs/main_config.yaml sim/output/
```

## Decision 1 — trajectory stability*
*Future work: Validate that stable ancestry trajectories can be inferred from simulated data*

## Decision 2 — null calibration (partial: matched-control neutral simulation)

`sim-ctrl-null` used `sim/validate_empirical_null.py --seed 251008` to exercise
`empirical_null` on **four independent neutral strata** (two chromosomes,
differentiated allele frequency, recombination/LD and observed time grids).
For each stratum, 48 neutral loci form a prespecified control pool and 32
independently simulated neutral loci are held out as candidate tests. Each
locus has three trajectories, four time bins and correlated shared neutral
drift across trajectories; counts are binomial at depths 70 or 130. QC-like
matching fields are drawn independently of outcomes. The matcher enforces
same chromosome, exact trajectory/time coverage, six numeric calipers and
at least 1 Mb distance, then selects 40 nearest neutral controls without
looking at scores. Every locus is analyzed with the same `trajectory_statistic`
call, including family concordance and GLS meta slope.

| Seed | Candidate loci | Controls/test | Mean empirical P | Uniform KS D | Fraction P ≤ 0.05 | Fraction P ≤ 0.10 |
|---|---:|---:|---:|---:|---:|---:|
| 251008 | 128 | 40 | 0.4825 | 0.1027 | 0.0625 | 0.1406 |

The deliberately broad automated bounds are mean P `[.43,.58]`, KS D `<.15`,
false-positive fraction at `.05 ≤ .10` and at `.10 ≤ .16`; all pass on this
seed. These are **single-seed descriptive checks**, not confidence intervals
or proof of uniformity. The discrete finite-control P grid (minimum `1/41`)
limits fine tail calibration. Stronger independent holdout simulations,
structure/replacement effects, locus-LD dependence, misspecified metadata,
trajectory-inference error, and actual ancient-DNA QC remain untested. Thus
Decision 2 is **not** a blanket go for biological selection claims. The
run's hashed inputs/outputs and package versions are in
`analysis/runs/20261008-0528-sim-ctrl-null-validation/manifest.json`.

## Decision 3 — power vs Akbari, FPR under migration (PRELIMINARY, underpowered)

Status: **in progress — directional only.** The run below used a deliberately fast
configuration (`n_loci = 12`, `runs = 6`, n ~= 66 tests per cell), so every number
carries binomial noise of roughly +/-0.05-0.08. It is a smoke-level result, not a
calibrated claim. A higher-n run is required before Decision 3 can be closed.

Method: for each scenario the FPR threshold is the 95th percentile of the score
distribution of *neutral* loci under the `neutral_drift` baseline, i.e. matched-null
calibration rather than an analytic Wald threshold. Three statistics are compared:
`trajectory_naive` (sign-concordance, no relatedness input), `trajectory_related`
(the same statistic with the relatedness matrix supplied so that the GLS
family-discount and `n_effective` are active), and `akbari` (the reconstructed
structure-corrected GLMM).

### Calibrated FPR (neutral / replacement)

| scenario | trajectory_naive | trajectory_related | akbari |
|---|---:|---:|---:|
| neutral_drift (baseline) | 0.061 | 0.061 | 0.083 |
| **replacement_only (baseline)** | **0.188** | **0.188** | **0.000** |

### Calibrated power (selection scenarios, baseline)

| scenario | trajectory_naive | trajectory_related | akbari |
|---|---:|---:|---:|
| one_population_selection | 0.500 | 0.667 | 0.167 |
| shared_selection | 1.000 | 1.000 | 1.000 |
| opposing_selection | 0.000 | 0.000 | 0.000 |
| selection_aligned_turnover | 1.000 | 1.000 | 1.000 |

### Answer to the central question

**[observation]** Under pure demographic replacement the naive replication statistic is
anti-conservative: calibrated FPR 0.188 versus 0.061 under neutral drift, while the
Akbari GLMM is conservative (0.000) in the same cell.

**[observation]** Supplying the relatedness matrix in this configuration changed the
replacement FPR by 0.000 (0.188 -> 0.188).

**[inference]** That null repair is expected, and is informative: the simulation
hands the statistic a *fixed* inter-trajectory relatedness of rho = 0.5, which yields
`n_effective` ~= 1.33 and therefore still credits two partly-independent
replications. In `replacement_only` one population literally becomes the other, so
its true relatedness is close to 1, `n_effective` should collapse toward 1, and a
single trajectory carries no concordance evidence. The load-bearing component of the
method is therefore **data-estimated** trajectory relatedness, not the concordance
count. Hand-set relatedness reproduces exactly the anti-conservative failure mode.

**[observation]** Where the GLMM loses power, the trajectory statistic recovers it:
`one_population_selection` calibrated power 0.667 (related) versus 0.167 (Akbari).

**[observation]** `opposing_selection` has power 0.000 for every statistic. This is the
expected design check: selection in opposing directions is not concordant, so a
replication statistic should be blind to it.

### Limits

Single fast configuration; rho supplied by hand rather than estimated; six scenarios
only; no linkage or LD structure among tested loci; no misspecification of the
trajectory-inference step beyond the `inference_error` factor; the FPR calibration
threshold is a finite-sample 95th percentile. Decision 3 is **not** closed.

Artifacts: `sim/sweeps/results/fpr_power.csv`, `sim/sweeps/results/summary.json`
(wall 62.6 s; seeds and thresholds recorded in the summary).

## Decision 4 — AMY and canonical loci*
*Future work: Validate method on known positive control loci*

## Third-pass RANKING validation — 2026-10-09 (append only)

[observation] Pre-run ranking metrics were frozen in docs/prioritization-prereg.md. Detailed completed results: docs/prioritization-results.md. SLC24A5 ranks9161/10641 and TYR1120 (outside decile); six/eight control regions reach the decile. Top100 control-proxy precision4%, region recall3/8; top10 proxy20%, recall2/8. Unknown selected loci are NOT validated negatives; these proxies are not actual modern-association precision.

[observation] At the observed top100 cutoff, separately simulated drift+jump loci exceed it32/300=10.67%, pure replacement18/300=6%, Brownian drift0/300 (95% upper1.22%). These are threshold tail costs, not FDR gates. Shape scores exact ramps1 and steps0, but noisy pure replacement median.267 versus drift.258 shows weak biological separation. Thirty node-stratified individual bootstraps retain43.7% of original top100 on average, Jaccard.281.

[inference] Positive-control enrichment supports exploratory candidate ranking; weak shape discrimination and moderate stability limit mechanistic claims. Selection elsewhere/earlier and demographic carriage remain valid reasons for follow-up, not exclusion. No offline modern trait enrichment was available/computed.
