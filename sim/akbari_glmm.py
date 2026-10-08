"""Akbari-style binomial time/GRM GLMM for small simulated data sets.

This is an independent, transparent PQL approximation to the model described in
Akbari et al. (2026) and PQLseqPy, NOT a numerical reproduction of its optimized
AI-PQL implementation. See docs/F-reproducible-code.md for assumptions and limits.
No randomness is used during fitting.
"""

from __future__ import annotations

from typing import Mapping

import numpy as np
import pandas as pd
from scipy import linalg, stats
from scipy.optimize import minimize
from scipy.special import expit


class AkbariGLMM:
    """Fit a per-locus binomial logit mixed model.

    y_i ~ Binomial(n_i, expit(beta_0 + beta_time * time_i + u_i)),
    u ~ N(0, tau1 * K + tau2 * I).  K is a *sample* GRM calculated from
    background markers excluding the tested locus's entire chromosome (LOCO).
    Both variance components are re-estimated for each locus by PQL/REML.

    Inputs are in sample order; each row is an observation at ONE locus, not a
    marker or a chromosome. Time units are supplied by the caller (e.g. years BP).
    Positive beta_time means increasing frequency toward *older* years BP.
    ``fit`` accepts either an externally constructed K or background genotypes;
    when target_chromosome is given with genotypes, marker_chromosomes is required.
    A supplied K must already exclude the target chromosome. No posterior
    probability is inferred: ``z_score`` is the two-sided Wald statistic.
    """

    def __init__(self, time_col="time_bp", allele_count_col="alt_count",
                 total_count_col="total_count"):
        self.time_col = time_col
        self.allele_count_col = allele_count_col
        self.total_count_col = total_count_col
        self.beta_time = None
        self.beta_intercept = None
        self.tau1 = None
        self.tau2 = None
        self.converged = False
        self.K = None

    @staticmethod
    def compute_grm(genotypes, method="vanraden", marker_chromosomes=None,
                    target_chromosome=None):
        """Return an n_samples x n_samples GRM from 0/1/2 background genotypes.

        Columns are centered on their own allele frequencies (NOT row means).
        Missing genotypes (NaN) are mean-imputed; invariant/zero-variance markers
        are discarded. ``vanraden`` uses sum_j 2p_j(1-p_j) as denominator;
        ``standardized`` divides each centered column by sqrt(2p_j(1-p_j))
        before averaging. No jitter is added. Chromosome labels are per COLUMN.
        """
        g = np.asarray(genotypes, dtype=float)
        if g.ndim != 2 or g.shape[0] < 4 or not np.all(np.isfinite(g) | np.isnan(g)):
            raise ValueError("genotypes must be a 2D finite/NaN sample-by-marker matrix")
        if np.any((g[~np.isnan(g)] < 0) | (g[~np.isnan(g)] > 2)):
            raise ValueError("genotypes must be in [0, 2]")
        if method not in ("vanraden", "standardized"):
            raise ValueError("method must be 'vanraden' or 'standardized'")
        if target_chromosome is not None and marker_chromosomes is None:
            raise ValueError("LOCO requires a chromosome label for every marker")
        if marker_chromosomes is not None:
            chrom = np.asarray(marker_chromosomes)
            if chrom.ndim != 1 or len(chrom) != g.shape[1]:
                raise ValueError("marker_chromosomes must have one label per column")
            if target_chromosome is not None:
                if not np.any(chrom == target_chromosome):
                    raise ValueError("target chromosome has no markers to exclude")
                g = g[:, chrom != target_chromosome]
        if g.shape[1] == 0:
            raise ValueError("no background markers remain")
        valid = np.any(~np.isnan(g), axis=0)
        g = g[:, valid]
        if not g.shape[1]:
            raise ValueError("no observed background markers remain")
        means = np.nanmean(g, axis=0)
        p = means / 2
        var = 2 * p * (1 - p)
        keep = var > 1e-10
        if not np.any(keep):
            raise ValueError("no polymorphic background markers remain")
        z = np.where(np.isnan(g[:, keep]), means[keep], g[:, keep]) - means[keep]
        if method == "vanraden":
            return (z @ z.T) / var[keep].sum()
        z /= np.sqrt(var[keep])
        return (z @ z.T) / z.shape[1]

    def _inputs(self, data, K):
        cols = (self.time_col, self.allele_count_col, self.total_count_col)
        if not all(c in data.columns for c in cols):
            raise ValueError(f"data requires columns {cols}")
        a = data.loc[:, cols].to_numpy(dtype=float)
        if not np.isfinite(a).all() or len(a) < 4:
            raise ValueError("at least four complete observations required")
        t, y, n = a.T
        if (np.any(n <= 0) or np.any(y < 0) or np.any(y > n)
                or np.any(n != np.floor(n)) or np.any(y != np.floor(y))):
            raise ValueError("counts must be integer with 0 <= alt <= total and total > 0")
        scale = np.std(t)
        if scale <= 0 or np.unique(t).size < 3:
            raise ValueError("time must have at least three distinct values")
        k = np.asarray(K, dtype=float)
        if k.shape != (len(t), len(t)) or not np.isfinite(k).all():
            raise ValueError("K must be a finite sample-by-sample matrix")
        if not np.allclose(k, k.T, atol=1e-8) or np.linalg.eigvalsh(k).min() < -1e-7:
            raise ValueError("K must be symmetric positive semidefinite")
        x = np.column_stack((np.ones(len(t)), (t - t.mean()) / scale))
        return x, y, n, (t.mean(), scale), (k + k.T) / 2

    @staticmethod
    def _gls(z, x, h):
        """REML ingredients for one PQL working Gaussian response."""
        c = linalg.cho_factor(h, check_finite=False)
        vi_x = linalg.cho_solve(c, x, check_finite=False)
        vi_z = linalg.cho_solve(c, z, check_finite=False)
        info = x.T @ vi_x
        ci = linalg.cho_factor(info, check_finite=False)
        cov = linalg.cho_solve(ci, np.eye(2), check_finite=False)
        beta = cov @ (x.T @ vi_z)
        residual = z - x @ beta
        py = linalg.cho_solve(c, residual, check_finite=False)
        reml = (2 * np.log(np.diag(c[0])).sum()
                + 2 * np.log(np.diag(ci[0])).sum() + residual @ py)
        return reml, beta, cov, py

    def fit(self, data: pd.DataFrame, K=None, genotypes=None,
            marker_chromosomes=None, target_chromosome=None, fixed_tau=None,
            max_iter=40, tol=1e-5):
        """Fit one locus and return beta_time, SE, Wald Z/p, tau1/tau2, convergence.

        ``fixed_tau=(tau1,tau2)`` is useful for a pre-estimated null covariance
        or the zero-variance GLM limit. Without it, nonnegative components are
        estimated by restricted working-Gaussian likelihood at each PQL step.
        A failed optimization raises RuntimeError, rather than returning a
        deceptively significant statistic. No rows are held out for LOCO.
        """
        if (K is None) == (genotypes is None):
            raise ValueError("provide exactly one of K and background genotypes")
        if K is None:
            K = self.compute_grm(genotypes, marker_chromosomes=marker_chromosomes,
                                 target_chromosome=target_chromosome)
        elif marker_chromosomes is not None:
            raise ValueError("marker_chromosomes only applies to background genotypes")
        x, y, n, (time_mean, time_scale), k = self._inputs(data, K)
        self.K = k
        if max_iter < 1 or tol <= 0:
            raise ValueError("max_iter and tol must be positive")
        if fixed_tau is not None:
            tau = np.asarray(fixed_tau, dtype=float)
            if tau.shape != (2,) or not np.isfinite(tau).all() or np.any(tau < 0):
                raise ValueError("fixed_tau must contain two nonnegative finite variances")
        else:
            tau = np.array([0.02, 0.02])
        # A pooled-count intercept starts GLM/PQL safely even for sparse alleles.
        beta = np.array([np.log((y.sum() + 0.5) / (n.sum() - y.sum() + 0.5)), 0.0])
        u = np.zeros(len(y))
        eye = np.eye(len(y))
        self.converged = False
        for iteration in range(1, max_iter + 1):
            eta = np.clip(x @ beta + u, -20, 20)
            p = expit(eta)
            w = np.maximum(n * p * (1 - p), 1e-8)
            z = eta + (y - n * p) / w
            r = 1 / w

            if fixed_tau is None:
                # Scale components to roughly unit diagonal K for optimizer.
                def objective(v):
                    try:
                        h = v[0] * k + np.diag(r + v[1])
                        return self._gls(z, x, h)[0]
                    except linalg.LinAlgError:
                        return np.inf

                opt = minimize(objective, tau, method="L-BFGS-B",
                               bounds=((0, None), (0, None)),
                               options={"maxiter": 70, "ftol": 1e-9})
                if not opt.success and not np.isfinite(opt.fun):
                    raise RuntimeError(f"variance optimization failed: {opt.message}")
                new_tau = np.maximum(opt.x, 0)
            else:
                new_tau = tau
            h = new_tau[0] * k + np.diag(r + new_tau[1])
            _, new_beta, cov, py = self._gls(z, x, h)
            new_u = (new_tau[0] * k + new_tau[1] * eye) @ py
            change = max(np.max(np.abs(new_beta - beta)),
                         np.max(np.abs(new_tau - tau)),
                         np.max(np.abs(new_u - u)))
            beta, tau, u = new_beta, new_tau, new_u
            if change < tol:
                self.converged = True
                break
        if not self.converged:
            raise RuntimeError(f"PQL did not converge in {max_iter} iterations")
        # Covariance of beta from final working weights and final variance model.
        p = expit(np.clip(x @ beta + u, -20, 20))
        h = tau[0] * k + np.diag(1 / np.maximum(n * p * (1 - p), 1e-8) + tau[1])
        _, _, cov, _ = self._gls(x @ beta, x, h)
        self.beta_time = float(beta[1] / time_scale)
        self.beta_intercept = float(beta[0] - beta[1] * time_mean / time_scale)
        self.tau1, self.tau2 = map(float, tau)
        se = float(np.sqrt(cov[1, 1]) / time_scale)
        z_score = self.beta_time / se
        return {"beta_time": self.beta_time, "beta_intercept": self.beta_intercept,
                "tau1": self.tau1, "tau2": self.tau2,
                "se_beta_time": se, "z_score": z_score,
                "p_value": float(2 * stats.norm.sf(abs(z_score))),
                "converged": True, "n_iterations": iteration}

    def leave_one_chromosome_out(self, loci: Mapping[str, pd.DataFrame],
                                 locus_chromosomes: Mapping[str, object],
                                 genotypes, marker_chromosomes):
        """Test each named locus using a background GRM excluding its chromosome.

        All locus data frames and genotype rows must refer to the SAME samples
        in the SAME order; locus_chromosomes maps each locus ID to a marker
        chromosome label. GRMs are cached once per target chromosome.
        Returns a dict of locus ID -> ``fit`` result. No test individuals are
        excluded. For genome-scale work, use precomputed LOCO K matrices.
        """
        if set(loci) != set(locus_chromosomes):
            raise ValueError("each locus needs exactly one chromosome label")
        cache = {}
        results = {}
        for name, data in loci.items():
            chrom = locus_chromosomes[name]
            if chrom not in cache:
                cache[chrom] = self.compute_grm(
                    genotypes, marker_chromosomes=marker_chromosomes,
                    target_chromosome=chrom)
            results[name] = self.fit(data, K=cache[chrom])
        return results


def simulate_data(n_individuals=80, n_timepoints=8, selection_strength=0.0,
                  population_structure=True, random_state=None, n_markers=120):
    """Seeded one-locus simulation; return (row-wise allele counts, background G).

    ``selection_strength`` is the log-odds change per SD of time BP, not a
    Wright-Fisher selection coefficient. Random effects are sampled from
    0.12*K + 0.04*I, with K from the returned (unlinked) background markers.
    Counts represent independent pooled allele trials conditional on the latent
    probability. This is a model-based check, not a demographic drift simulator.
    """
    if random_state is None:
        raise ValueError("simulation requires explicit random_state seed")
    if n_individuals < 4 or n_timepoints < 3 or n_markers < 2:
        raise ValueError("simulation requires >=4 samples, >=3 times, >=2 markers")
    rng = np.random.default_rng(random_state)
    t = rng.choice(np.linspace(0, 10000, n_timepoints), n_individuals)
    scaled = (t - t.mean()) / t.std()
    ancestry = rng.binomial(1, expit(0.7 * scaled)) if population_structure else np.zeros(n_individuals)
    base = rng.uniform(0.2, 0.8, n_markers)
    shift = rng.uniform(-0.15, 0.15, n_markers) if population_structure else np.zeros(n_markers)
    g = rng.binomial(2, np.clip(base[None, :] + ancestry[:, None] * shift, 0.01, 0.99))
    k = AkbariGLMM.compute_grm(g)
    u = linalg.cholesky(0.12 * k + 0.04 * np.eye(n_individuals), lower=True) @ rng.standard_normal(n_individuals)
    depth = rng.integers(8, 17, size=n_individuals)
    y = rng.binomial(depth, expit(-0.8 + selection_strength * scaled + u))
    data = pd.DataFrame({"time_bp": t, "alt_count": y, "total_count": depth})
    return data, g
