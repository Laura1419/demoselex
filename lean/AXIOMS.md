# AXIOMS Documentation

This file documents all uses of `sorry`, `axiom`, and `native_decide` in the Lean codebase according to the honest accounting policy.

## Policy

Every use of the following constructs must be documented here:
- `sorry` - indicates incomplete proofs needing completion
- `axiom` - foundational assumptions requiring justification  
- `native_decide` - computationally intensive decidability proofs

## Current Axioms and Gaps

### Module: Core.lean
**Status**: No axioms or gaps found

### Module: Trajectory.lean
**Status**: No axioms or gaps found

### Module: FrequencyPath.lean
**Gaps identified**:
- `concordance_score` function uses `pure 0.5` placeholder
  - **Location**: Line 22
  - **Justification**: Need to define formal concordance criteria
  - **Priority**: Medium - affects core analysis

### Module: Concordance.lean  
**Gaps identified**:
- `check_replication` uses placeholder threshold of 0.7
  - **Location**: Line 19
  - **Justification**: Statistical threshold needs empirical validation
  - **Priority**: High - affects significance determination
- `statistical_significance` returns fixed value 0.05
  - **Location**: Line 23
  - **Justification**: Needs proper statistical test implementation
  - **Priority**: High - critical for scientific validity

## Completion Plan

### Phase 1: Core Definitions ✓
- [x] Basic type definitions
- [x] Trajectory data structure
- [x] Direction enumeration

### Phase 2: Analysis Functions ⚠️
- [x] Direction calculation
- [x] Monotonicity checking
- [ ] **PENDING**: Concordance scoring algorithm
- [ ] **PENDING**: Statistical significance testing

### Phase 3: Validation
- [ ] Simulation-based validation
- [ ] Statistical calibration
- [ ] Empirical threshold determination

## Notes

- Current implementation focuses on structural correctness
- Statistical methods need integration with Mathlib's probability theory
- Thresholds should be determined through null simulations
- The 0.7 concordance threshold is provisional and requires validation