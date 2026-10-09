#!/usr/bin/env python3
"""Exact algebra checks for the preregistration; NOT an ancestry painter/estimator.

All examples are synthetic. No genotype bundle or AMY outcome is loaded here.
"""
from __future__ import annotations
import numpy as np
from scipy.stats import norm

AMY_EXCLUSION = (103915000, 104327000)  # inclusive GRCh37 positions


def genotype_probabilities(p_a: float, p_b: float) -> np.ndarray:
    """Unphased two-Bernoulli convolution, including a=b diploid states."""
    p = np.asarray([p_a, p_b], dtype=float)
    if not np.isfinite(p).all() or (p < 0).any() or (p > 1).any():
        raise ValueError("Allele frequencies must be finite and in [0,1]")
    return np.array([(1-p_a)*(1-p_b), p_a*(1-p_b)+(1-p_a)*p_b, p_a*p_b])


def pair_marginal_distribution(frequencies, pair_weights) -> np.ndarray:
    """Marginalize an UNORDERED pair distribution inside the likelihood.

    pair_weights maps (a,b), a<=b, to calibrated pair probability; this must
    not be reconstructed by multiplying marginal ancestry dosages.
    """
    p = np.asarray(frequencies, dtype=float)
    if p.ndim != 1 or not len(p):
        raise ValueError("Need a nonempty frequency vector")
    if not pair_weights or not np.isclose(sum(pair_weights.values()), 1., atol=1e-10):
        raise ValueError("Pair probabilities must sum to one")
    distribution = np.zeros(3)
    for (a, b), weight in pair_weights.items():
        if not 0 <= a <= b < len(p) or not np.isfinite(weight) or weight < 0:
            raise ValueError("Invalid unordered state pair/probability")
        distribution += weight * genotype_probabilities(p[a], p[b])
    return distribution


def shared_direction_statistic(contrasts) -> float:
    """Smallest common-direction contrast; nonpositive for discordance/zero."""
    delta = np.asarray(contrasts, dtype=float)
    if delta.ndim != 1 or len(delta) < 2 or not np.isfinite(delta).all():
        raise ValueError("Require finite contrasts from at least two populations")
    return float(max(delta.min(), (-delta).min()))


def spectral_effective_count(correlation) -> float:
    """Independent score dimensions, NOT a count of genealogical ancestors."""
    r = np.asarray(correlation, dtype=float)
    if r.ndim != 2 or r.shape[0] != r.shape[1] or not len(r):
        raise ValueError("Need a nonempty square correlation matrix")
    if not np.isfinite(r).all() or not np.allclose(r, r.T, atol=1e-10):
        raise ValueError("Correlation must be finite and symmetric")
    if not np.allclose(np.diag(r), 1., atol=1e-10) or np.linalg.eigvalsh(r).min() < -1e-10:
        raise ValueError("Correlation must be positive semidefinite with unit diagonal")
    return float(np.trace(r)**2 / np.sum(r*r))


def ideal_frequency_floor(n_early: int, n_late: int, p: float = .5) -> float:
    """Optimistic independent-diploid two-bin MDE; not achieved ancient power."""
    if n_early < 1 or n_late < 1 or not np.isfinite(p) or not 0 < p < 1:
        raise ValueError("Need positive bin counts and a polymorphic frequency")
    return float((norm.ppf(.975)+norm.ppf(.8)) *
                 np.sqrt(p*(1-p)/(2*n_early)+p*(1-p)/(2*n_late)))


def outside_amy_exclusion(chromosomes, positions) -> np.ndarray:
    """Coordinate exclusion only; this does NOT certify imputation independence."""
    chrom, pos = np.asarray(chromosomes), np.asarray(positions)
    if chrom.shape != pos.shape:
        raise ValueError("Chromosome/position shapes differ")
    lo, hi = AMY_EXCLUSION
    return ~((chrom == 1) & (pos >= lo) & (pos <= hi))
