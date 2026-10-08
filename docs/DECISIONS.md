# Decision Points

Critical go/no-go decisions for the trajectory-replication research project.

| Decision | Question | If No | Responsible Tasks | Evidence Required | Status |
|----------|----------|-------|-------------------|-------------------|--------|
| **1** | Can stable ancestry trajectories actually be inferred? | reformulate | `traj-validate` | traj-validation.md shows stable structure, traj-quality-report | PENDING |
| **2** | Under neutral demographic simulations, is the replication statistic calibrated? | fix before touching biological loci | `sim-sweeps`, `sim-ctrl-null` | neutral FPR ≤ 0.05, calibration plots | PENDING |
| **3** | Does the method recover simulated repeated selection better than Akbari without greatly increasing false positives? | the method has not solved the problem | `sim-sweeps` | power comparison, ROC curves, effect-size bias | PENDING |
| **4** | Do AMY and canonical loci behave sensibly? | investigate before scaling | `locus-amy`, `locus-canonical` | AMY trajectory patterns, known locus replication | PENDING |
| **5** | Does the Akbari reanalysis reveal a meaningful low-Akbari/high-replication class? | conclude the hypothesised over-correction appears limited under this test | `locus-akbari` | class exists and replicates, enrichment statistics | PENDING |

## Decision Tracking

Each decision will be recorded here with:
- **Status**: PENDING / IN_PROGRESS / GO / NO_GO
- **Date**: When the decision was made
- **Evidence**: Link(s) to the supporting analysis/documentation
- **Rationale**: Brief explanation of the reasoning

This file should be updated as each decision point is reached and resolved.