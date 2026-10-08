# F. Reproducible code

> Owner tasks: all. Status: in progress (Lean bootstrap).

## Lean formalization (lean-setup)

The canonical Lake package lives in `lean/` (`lean/lakefile.lean`, `lean/lean-toolchain`, `lean/lake-manifest.json`). A repository-root Lake entry point also builds the same `lean/TrajectorySelection.lean` and its `TrajectorySelection/{Core,Trajectory,FrequencyPath,Concordance}.lean` imports. Both entry points pin Lean 4.34.1 and Mathlib v4.34.1; the checked-in manifests lock Mathlib and transitive package revisions. Run `cd lean && lake exe cache get && lake build` or, from the repository root, `lake exe cache get && lake build`. Cache downloads are optional for correctness; building Mathlib from source uses more time and memory.

`lean/AXIOMS.md` is the mandatory ledger for every project `sorry`, `axiom`, and `native_decide` and flags inherited scientific placeholders separately. A successful build validates import/elaboration only, not scientific calibration. Definitions and proof attempts are owned by downstream Lean tasks; do not interpret placeholder scores as selection evidence.

The general run-manifest harness is `scripts/provenance.py` (see `docs/PROVENANCE.md`); a Lean build records seed `0` as a no-randomness convention and captures toolchain and Mathlib revision with the inputs in `analysis/runs/`.

Further pipeline script inventory (downloads/QC/trajectories/loci/simulations/figures/reanalysis) is pending its owner tasks.

## Trajectory replication statistic (`sim-stat`)

`sim/trajectory_stat.py` exposes `trajectory_statistic(data, *, family=None,
relatedness=None) -> TrajectoryResult`. Supply a pandas long table with
`locus,trajectory,time,alt_count,total_count`; see Section D for the equations,
family collapse, assumptions and interpretation. The returned `.slopes` is a
tidy per-trajectory table; `.summary` has the pooled Model 0 slope, raw
concordance, relatedness-discounted Model 2 GLS slope/SE, effective family
count and heterogeneity diagnostics. No RNG occurs in fitting. A minimal call:

```python
from trajectory_stat import trajectory_statistic
result = trajectory_statistic(counts, family={"a": "a", "b": "b"},
                              relatedness=family_correlation)
print(result.slopes[["locus", "trajectory", "slope", "se", "direction", "n", "span"]])
print(result.summary[["locus", "concordance", "n_effective", "meta_slope", "q"]])
```

From the repository root, run `python3 -m pytest -q sim/test_trajectory_stat.py`
and `python3 -m sim.validate_trajectory_stat --seed 31415`. The latter records
the seeded synthetic counts, tidy slopes, summary, code hashes and software
versions through `scripts/provenance.py` under `analysis/runs/`.
The fixed-seed designed test checks strong trends of both signs; analytic
binomial SE, duplication invariance, relatedness discount, zero/tied signs,
heterogeneity and a turnover-only pooled trend are separate regressions.
Assumed binomial independence and family-error correlation must be evaluated
under the matched null in downstream simulation; diagnostic p-values are not
yet calibrated genome-wide significance.

## Akbari-style GLMM comparator (`sim-glmm`)

`sim/akbari_glmm.py` implements an independent small-sample binomial PQL/REML
approximation of Akbari et al. (2026)'s PQLseqPy model, not the production scan
or an estimate of a Wright–Fisher selection coefficient. Per locus,
`y_i ~ Binomial(n_i, expit(beta_0 + beta_time*t_i + u_i))` and
`u ~ N(0, tau1*K_LOCO + tau2*I)`. Working response
`z = eta + (y-n*p)/(n*p*(1-p))`; working covariance
`H = diag(1/(n*p*(1-p))) + tau1*K_LOCO + tau2*I`.
Nonnegative variance components are optimized by restricted working-Gaussian
likelihood **separately for each locus**, then the fixed effects and BLUP are
updated. The covariance `(X.T @ H^-1 @ X)^-1` supplies a two-sided Wald Z/p.
Time is centered/scaled internally; returned `beta_time` and SE are per unit
of input time. For years BP, a positive slope denotes higher frequency in
*older* samples, not an increase toward the present. No posterior selection
probability is claimed.

Minimal seeded example from the repository root (requires numpy, pandas, scipy):

```bash
PYTHONPATH=sim python3 - <<'PY'
from akbari_glmm import AkbariGLMM, simulate_data
import numpy as np
counts, background = simulate_data(random_state=42, selection_strength=1.2)
chrom = np.repeat(["1", "2", "3", "4"], 30)  # one label per background column
result = AkbariGLMM().fit(counts, genotypes=background,
                          marker_chromosomes=chrom, target_chromosome="1")
print(result["z_score"], result["p_value"], result["converged"])
PY
python3 -m pytest sim/test_akbari_glmm.py -q
```

For many target loci, use `leave_one_chromosome_out(loci, locus_chromosomes,
background_genotypes, marker_chromosomes)` where `loci` maps locus IDs to
sample-ordered count frames and locus_chromosomes maps IDs to labels. The
method caches **one K per excluded marker chromosome**, retains all sample
rows, and fits each locus separately. Alternatively pass a precomputed
sample-by-sample LOCO `K` to `fit`; the caller then owns exclusion. The GRM
centers each marker by its allele frequency, mean-imputes missing calls and
omits invariant markers. A chromosome label belongs to a marker column, not
an individual. No genotype dosage/imputation uncertainty is modeled.

**Validation and external check.** `sim/test_akbari_glmm.py` tests a 30-seed
structured null (no skipped fits, mean |Z| bound, SD and false-positive bound),
a 10-seed strong single-locus time trend (each exceeds |Z|=5.45), the GLM
Fisher-information limit, and LOCO exclusion. On the seeded 48-sample fixture
with fixed `(tau1,tau2)=(0.12,0.04)`, our Z values -0.768104 (null) and 8.142206
(selected) agree within 0.003 with the **released PQLseqPy** implementation
at GitHub `mokar2001/PQLseqPy` commit `0934722dfcc4bae1cab201ce6db8cd20e08a22dd`
(`GLMM(X,Y,K,fixed_tau=(.12,.04)).fit()`, `X=[1, standardized_time]`,
`Y=[alt, total-alt]`). This reference comparison was executed with isolated
statsmodels 0.14.5/patsy 1.0.1; regression fixtures do not depend on them.
The full 9.7M-variant published summary cannot be compared locus-by-locus
without matching sample-count matrices, GRMs and QC; agreement on two
simulated fits is **not** empirical validation of Akbari's reported loci or
paper-figure FPR. Small-n PQL Wald tails can be biased, particularly at
rare alleles or time/ancestry collinearity. Downstream `sim-sweeps` must assess
calibration, power and replacement nulls in the demographic simulator before
any biological interpretation.
