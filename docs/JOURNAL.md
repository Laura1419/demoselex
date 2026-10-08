# Research journal

Append-only. Newest entries first. Format:

```
## YYYY-MM-DD HH:MMZ — <task-id> — <one-line>
- what was done
- what was found
- blockers / next
```

## 2026-10-08 04:05Z — recon-akbari — Comprehensive Akbari methodology reconstruction completed
- Read Akbari et al. 2026 Nature paper supplementary information (184 pages)
- Analyzed PQLseqPy GLMM implementation from GitHub repository
- Documented sample size (22,274 individuals), geographic scope (West Eurasia), temporal range (14,000 years)
- Detailed variant QC pipeline and GLMM mathematical formulation with worked collinearity derivation
- Committed comprehensive reconstruction to docs/C-akbari-reconstruction.md with seed tracking note
- Found that provenance.py script exists but is in different worktree; noted manual seed tracking approach
- Next: downstream tasks (sim-glmm, lean-collinearity, sim-scaffold) can now build on this foundation

## 2026-10-08 05:00Z — sim-scaffold — Repaired and validated dated mosaic simulator
- Replaced the earlier final-only, non-donor migration prototype with sampled diploid mosaic Wright–Fisher trajectories, actual donor migration/turnover, and signed selection.
- All six fixed-seed configs now emit aligned dated genotype/sample matrices, per-trajectory locus truth, and `scripts/provenance.py` manifests; end-to-end tests pass. SLiM was unavailable.
- Pure replacement has zero selection truth despite a locus-frequency rise; opposing selection has negative truth in pop2. These are scenario controls, not power estimates.
- Next: sim-sweeps should vary seeds, migration and turnover to measure calibration, power and model misspecification.
