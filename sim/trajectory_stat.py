"""Within-trajectory frequency trends and relatedness-discounted locus summaries.

The input is a long table of independent binomial allele-count observations with
columns ``locus, trajectory, time, alt_count, total_count``. ``time`` is numeric:
a positive slope means an increase with *increasing numeric time* (the opposite
of an increase toward the present if time is years BP). Observations sharing a
locus, trajectory and time are pooled before fitting. The estimator is the
sample-size-weighted frequency slope from docs/D-method.md, not a logit slope.

The SE assumes independent binomial sampling of each time bin, with a Jeffreys
plug-in (alt+1/2)/(total+1) to avoid a spurious zero SE for boundary counts.
It does not model drift, date uncertainty, batch effects or overdispersion.
Relatedness is a *prespecified correlation of family-level slope errors*, not
an ancestry matrix automatically estimated from the tested locus. The Wald and
heterogeneity p-values are working-Gaussian diagnostics; calibrate genomic
significance with matched controls and neutral demographic simulations.
No randomness occurs during estimation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd
from scipy import linalg, stats


INPUT_COLUMNS = ("locus", "trajectory", "time", "alt_count", "total_count")
SLOPE_COLUMNS = ("locus", "trajectory", "slope", "se", "direction", "n", "span")


@dataclass(frozen=True)
class TrajectoryResult:
    """Tidy slopes (one row/trajectory/locus) and locus-level GLS summaries."""

    slopes: pd.DataFrame
    summary: pd.DataFrame


def _concordance(values: np.ndarray) -> float:
    """Max(positive,negative)/informative; undefined for all-zero slopes."""
    positive = int(np.count_nonzero(values > 0))
    negative = int(np.count_nonzero(values < 0))
    return max(positive, negative) / (positive + negative) if positive + negative else float("nan")


def _validate(data: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(data, pd.DataFrame) or not set(INPUT_COLUMNS).issubset(data.columns):
        raise ValueError(f"data requires columns {INPUT_COLUMNS}")
    if data.empty or data.loc[:, INPUT_COLUMNS].isna().any().any():
        raise ValueError("data must contain nonmissing observations")
    rows = data.loc[:, INPUT_COLUMNS].copy()
    for name in ("locus", "trajectory"):
        if not rows[name].map(lambda x: isinstance(x, str) and bool(x.strip())).all():
            raise ValueError(f"{name} labels must be nonempty strings")
    try:
        values = rows.loc[:, ("time", "alt_count", "total_count")].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("time and counts must be numeric") from exc
    if not np.isfinite(values).all():
        raise ValueError("time and counts must be finite")
    time, alt, total = values.T
    if (np.any(total <= 0) or np.any(alt < 0) or np.any(alt > total)
            or np.any(alt != np.floor(alt)) or np.any(total != np.floor(total))):
        raise ValueError("counts must be integer with 0 <= alt <= total and total > 0")
    rows.loc[:, ("time", "alt_count", "total_count")] = values
    return rows


def _weighted_slope(group: pd.DataFrame, label: str) -> tuple[float, float, int, float]:
    bins = group.groupby("time", as_index=False, sort=True)[["alt_count", "total_count"]].sum()
    t = bins.time.to_numpy(dtype=float)
    y = bins.alt_count.to_numpy(dtype=float)
    w = bins.total_count.to_numpy(dtype=float)
    if len(t) < 2:
        raise ValueError(f"{label} requires at least two distinct times")
    centered = t - np.average(t, weights=w)
    sxx = float(np.dot(w, centered * centered))
    if not np.isfinite(sxx) or sxx <= 0:
        raise ValueError(f"{label} has undefined weighted time variation")
    # Since sum(w*(t-tbar))=0, the intercept term cancels in the numerator.
    slope = float(np.dot(centered, y) / sxx)
    # Floating cancellation on a mathematically constant frequency path can
    # yield a spurious ~1e-16 sign. This is only a roundoff guard, not a
    # statistical significance threshold for small but real trends.
    roundoff = (32 * np.finfo(float).eps * max(float(np.max(y / w)), 1.0)
                / float(t.max() - t.min()))
    if abs(slope) <= roundoff:
        slope = 0.0
    # Var(sum_i a_i * (y_i/w_i)) = sum_i a_i^2 p_i(1-p_i)/w_i,
    # where a_i=w_i*(t_i-tbar)/Sxx and counts in distinct bins are independent.
    p = (y + 0.5) / (w + 1.0)
    se = float(np.sqrt(np.dot(w * centered * centered, p * (1 - p))) / sxx)
    if not np.isfinite(se) or se <= 0:
        raise ValueError(f"{label} has undefined binomial SE")
    return slope, se, len(t), float(t.max() - t.min())


def _fit_slopes(data: pd.DataFrame) -> pd.DataFrame:
    records = []
    # Grouping/sorting also makes results independent of the order of input rows.
    for (locus, trajectory), group in data.groupby(["locus", "trajectory"], sort=True):
        slope, se, n, span = _weighted_slope(group, f"{locus}/{trajectory}")
        records.append((locus, trajectory, slope, se, int(np.sign(slope)), n, span))
    return pd.DataFrame.from_records(records, columns=SLOPE_COLUMNS)


def _correlation(relatedness: pd.DataFrame | None, families: list[str]) -> np.ndarray:
    if relatedness is None:
        return np.eye(len(families))  # Explicit working independence assumption.
    if (not isinstance(relatedness, pd.DataFrame) or not relatedness.index.is_unique
            or not relatedness.columns.is_unique or
            not set(families).issubset(relatedness.index) or
            not set(families).issubset(relatedness.columns)):
        raise ValueError("relatedness must be a labeled matrix covering every family")
    try:
        r = relatedness.loc[families, families].to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("relatedness entries must be numeric") from exc
    if (not np.isfinite(r).all() or not np.allclose(r, r.T, atol=1e-10)
            or not np.allclose(np.diag(r), 1, atol=1e-10)
            or np.any(r < -1e-10) or np.any(r > 1 + 1e-10)):
        raise ValueError("relatedness must be symmetric with diagonal 1 and entries in [0,1]")
    # An exact duplicate belongs to the same family; do not regularize a
    # singular matrix into apparently independent evidence.
    try:
        linalg.cholesky(r, lower=True)
    except linalg.LinAlgError as exc:
        raise ValueError("relatedness must be positive definite; merge identical families") from exc
    return r


def trajectory_statistic(
    data: pd.DataFrame,
    *,
    family: Mapping[str, str] | None = None,
    relatedness: pd.DataFrame | None = None,
) -> TrajectoryResult:
    """Fit weighted slopes and per-locus correlated fixed-effect GLS meta-analysis.

    ``family`` maps each trajectory to a *prespecified* relatedness family.
    If omitted, each trajectory is its own family. Copies in one family are
    combined by an equal-weight mean and assigned SE=max(member SEs, sample
    SD of member slopes); thus copying an identical path does not increase
    precision or N_eff. This conservative collapse is a sensitivity model,
    not a learned clustering algorithm or a likelihood for within-family data.
    ``relatedness`` is a pandas correlation DataFrame indexed by the family
    labels, with PSD strengthened to positive definite for identifiable GLS.
    If omitted, distinct families are assumed independent. It must be set
    from background evidence rather than estimated at the tested locus.

    ``slopes`` has the required locus/trajectory/slope/se/direction/n/span
    columns (n is distinct time bins; span is in supplied time units). Zero
    slopes have direction 0 and are omitted from concordance denominators.
    ``summary.pooled_slope`` is the unstratified Model 0 count trend, for
    comparison only: turnover can create this trend without within-lineage
    change. ``summary.concordance`` follows the Lean trajectory sign definition;
    ``family_concordance`` discounts copies by using collapsed family signs.
    ``n_effective`` is m²/sum(R) on distinct families, not a fractional
    binomial trial count. ``meta_slope`` and ``meta_se`` use V=diag(SE) R
    diag(SE), and Q is the corresponding residual quadratic form (df=m-1).
    An all-zero concordance or single-family heterogeneity is undefined (NaN).
    """
    data = _validate(data)
    slopes = _fit_slopes(data)
    trajectories = set(slopes.trajectory)
    if family is None:
        families = {trajectory: trajectory for trajectory in trajectories}
    else:
        if set(family) != trajectories or any(
            not isinstance(v, str) or not v.strip() for v in family.values()
        ):
            raise ValueError("family must map every trajectory to a nonempty family label")
        families = dict(family)
    slopes["family"] = slopes.trajectory.map(families)
    summaries = []
    for locus, group in slopes.groupby("locus", sort=True):
        pooled, pooled_se, _, _ = _weighted_slope(data.loc[data.locus == locus], locus)
        collapsed = []
        for label, members in group.groupby("family", sort=True):
            effects = members.slope.to_numpy(dtype=float)
            errors = members.se.to_numpy(dtype=float)
            spread = float(np.std(effects, ddof=1)) if len(effects) > 1 else 0.0
            collapsed.append((label, float(effects.mean()), max(float(errors.max()), spread)))
        labels = [row[0] for row in collapsed]
        beta = np.array([row[1] for row in collapsed])
        se = np.array([row[2] for row in collapsed])
        r = _correlation(relatedness, labels)
        v = se[:, None] * r * se[None, :]
        cho = linalg.cho_factor(v, check_finite=False)
        one = np.ones(len(beta))
        inv_one = linalg.cho_solve(cho, one, check_finite=False)
        info = float(one @ inv_one)
        meta = float((beta @ inv_one) / info)
        meta_se = float(np.sqrt(1 / info))
        residual = beta - meta
        q = float(residual @ linalg.cho_solve(cho, residual, check_finite=False))
        df = len(beta) - 1
        summaries.append({
            "locus": locus, "pooled_slope": pooled, "pooled_se": pooled_se,
            "n_trajectories": len(group), "n_families": len(beta),
            "n_effective": float(len(beta)**2 / r.sum()),
            "concordance": _concordance(group.slope.to_numpy(dtype=float)),
            "family_concordance": _concordance(beta),
            "meta_slope": meta, "meta_se": meta_se,
            "z_score": meta / meta_se,
            "p_value": float(2 * stats.norm.sf(abs(meta / meta_se))),
            "q": q if df else float("nan"), "q_df": df,
            "q_p_value": float(stats.chi2.sf(q, df)) if df else float("nan"),
            "i_squared": max(0.0, (q - df) / q) if df and q > 0 else float("nan"),
        })
    return TrajectoryResult(slopes=slopes, summary=pd.DataFrame(summaries))
