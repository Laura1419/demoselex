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

## Decision 2 — null calibration*
*Future work: Under neutral demographic simulations, calibrate the replication statistic*

## Decision 3 — power vs Akbari, FPR under migration*
*Future work: Compare method performance against Akbari baseline*

## Decision 4 — AMY and canonical loci*
*Future work: Validate method on known positive control loci*
