# Decision Points

Critical go/no-go decisions for the trajectory-replication research project.
Updated 2026-10-08 after the Tier-1 real-data reanalysis, the calibrated subset
simulation, and the Tier-2 data acquisition.

| Decision | Question | Status | Evidence |
|---|---|---|---|
| **1** | Can stable ancestry trajectories actually be inferred? | **GO (partly)** | `docs/traj-validation.md`, `docs/traj-v3-merge.md`. v1 failed (one 9,648-individual blob spanning 65,000 yr). v2 recovered 13 label-free paths and WHG-like / early-farmer-like / steppe-like histories emerged post-hoc from genotypes alone. Caveat: Avar-period over-split and a 4,594-individual under-resolved blob; honest replication count **N_eff ~ 3-5, use ~4**. |
| **2** | Under neutral demographic simulations, is the replication statistic calibrated? | **GO with a stated cost** | `docs/G-validation-report.md` Decision 2, `sim/validate_empirical_null.py`. Matched-control calibration holds under `neutral_drift` (FPR 0.05). |
| **3** | Does the method recover simulated repeated selection better than Akbari without greatly increasing false positives? | **CONDITIONAL GO - regime-specific** | `docs/G-validation-report.md` Decision 3, `docs/sim-subset-power.md`. The method does **not** beat the GLMM in general. It beats it decisively in one regime: selection in a **minority of trajectories** (1 of 6) -> calibrated power **0.95 vs 0.25** (strength 1.0). Costs: calibrated FPR under pure replacement **0.112 vs 0.045**; the max-over-subsets device is worthless uncalibrated (FPR 0.55-0.64) and its argmax is the maximum allowed size for 87% of loci. |
| **4** | Do AMY and canonical loci behave sensibly? | **PARTIAL** | `docs/H-amy-analysis.md`. AMY: 537 ancient carriers, 522 map to AADR with dates+coordinates, 288 across 10 structural haplotypes; the five AMY questions are blocked on trajectories + sequence-level and relatedness data. Canonical-locus check not yet run against the v3 units. |
| **5** | Does the Akbari reanalysis reveal a meaningful low-Akbari/high-replication class? | **NO (Tier-1); Tier-2 decisive test in flight** | `docs/I-quadrant-first-pass.md`, `docs/subset-regime-results.md`. Low-\|X\| loci show high replication at **4.27%** vs the 5% matched-null rate (no excess); high-\|X\| loci at **9.24%** (~2x). The pre-registered subset-regime test found **0/2,000 at BH-FDR 0.05** (min q 0.639), candidate enrichment **1.13x**. Tier-2 (same individuals/variants as the paper) is running. |

## Decision Tracking

### Decision 1 - rationale
Trajectory inference was iterated once (v1 -> v2) rather than accepted. v1's failure was
diagnosed (moderns glueing the graph; no within-trajectory continuity requirement) and
fixed with an ordered-path construction. v3 merged relatedness groups but *does not*
yield a defensible large replication count: there is no natural break in the relatedness
spectrum, and G0 still bundles 10,058 individuals of mixed period/region origin. The
load-bearing number is therefore N_eff ~ 4, not 13 and not 6.

### Decision 3 - rationale
This is the most important correction to the original framing. There is no general
power advantage. The advantage is confined to the *minority-selection* regime, is
quantified (0.95 vs 0.25 at f = 1/6, strength 1.0), and carries a measured FPR price.
Where selection is carried in by an expanding population, **both** methods are
conservative: Akbari because the time vector is collinear with the ancestry axis, ours
because within-trajectory change is flat for a carried allele. The two blind spots
overlap, which is the structural reason the Tier-1 pooled comparison shows no
over-correction class.

### Decision 5 - rationale
The first pass initially returned Q3 = 0 *by construction* (candidates were defined as
|X| >= 3.16 and the empirical null was computed only for candidates). After rebuilding a
leave-one-out matched null for the low-|X| pool, Q3 becomes well-posed at 635 loci -
which is 4.3% of the pool, i.e. **at or below** the 5% null rate. Counting loci at an
empirical-P threshold yields ~5% by construction; only a departure from that rate is
evidence, and there is none. Akbari-significant loci are ~2x *more* likely to also show
high replication, so the two statistics agree rather than dissociate.

### Open items
- Tier-2 scan (`analysis/tier2_scan.py`) on 5,788 Akbari-analysed individuals and the
  paper's own variant set - the decisive test of Decision 5.
- Canonical-locus check (LCT, SLC24A5/SLC45A2/TYR/OCA2, FADS1/2, immune, AMY) against
  the v3 units, to test the parallel-adaptation prior directly.
- A better-behaved subset statistic (pre-specified subset size or complement-flatness
  score) instead of the unconstrained max.
- `lean-review` audit: `lean/AXIOMS.md` still lists three provisional placeholder
  definitions that must not be cited.
