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

## Decision 3 — power vs Akbari, FPR under migration*
*Future work: Compare method performance against Akbari baseline*

## Decision 4 — AMY and canonical loci*
*Future work: Validate method on known positive control loci*
