"""Matched-locus empirical null for within-trajectory replication.

Supply *prespecified*, QC-derived locus metadata and allele-count observations to
``empirical_null``. Matching never uses a locus's trajectory statistic. All
candidate and control loci are fitted by the same ``trajectory_statistic`` call;
no analytic Wald p-value is interpreted as demographic calibration.

Time bins must be harmonized across loci. The exact set of (trajectory, time)
observations is required to match, in addition to the supplied temporal-coverage
fraction. Same-chromosome controls closer than ``exclusion_bp`` are removed to
reduce local LD leakage. Remaining controls must satisfy every absolute caliper;
nearest controls are chosen by a prespecified scaled L1 distance, with locus ID
as deterministic tie-break. Defaults are illustrative, not universal QC rules.
Failed support/calipers/control-count checks raise ValueError, never a fallback
unmatched null. The caller must exclude selected or linked loci from the
background pool and audit ancestry, batch effects, call quality, and genomic
independence before biological inference.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

try:  # Both ``python3 -m sim...`` and PYTHONPATH=sim imports.
    from .trajectory_stat import trajectory_statistic
except ImportError:
    from trajectory_stat import trajectory_statistic


FEATURES = ("allele_frequency", "recombination_rate", "local_ld",
            "imputation_quality", "missingness", "temporal_coverage")
# Recombination rate is in cM/Mb; other features are fractions/proportions.
DEFAULT_CALIPERS = {
    "allele_frequency": 0.05, "recombination_rate": 0.5,
    "local_ld": 0.1, "imputation_quality": 0.05,
    "missingness": 0.05, "temporal_coverage": 0.05,
}
METADATA_COLUMNS = ("locus", "chromosome", "position", *FEATURES)


@dataclass(frozen=True)
class EmpiricalNullResult:
    """Candidate ranks, per-match audit rows, and identical fitted statistics."""

    results: pd.DataFrame
    matches: pd.DataFrame
    statistics: pd.DataFrame


def _metadata(frame: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame) or not set(METADATA_COLUMNS) <= set(frame):
        raise ValueError(f"metadata requires columns {METADATA_COLUMNS}")
    if frame.empty or frame.loc[:, METADATA_COLUMNS].isna().any().any():
        raise ValueError("metadata must be nonempty and complete")
    table = frame.loc[:, METADATA_COLUMNS].copy()
    for label in ("locus", "chromosome"):
        if not table[label].map(lambda x: isinstance(x, str) and bool(x.strip())).all():
            raise ValueError(f"{label} must contain nonempty string labels")
    if table.locus.duplicated().any():
        raise ValueError("metadata requires one row per locus")
    for label in ("position", *FEATURES):
        try:
            table[label] = pd.to_numeric(table[label], errors="raise").astype(float)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{label} must be numeric") from exc
        if not np.isfinite(table[label].to_numpy()).all():
            raise ValueError(f"{label} must be finite")
    if ((table.position < 0) | (table.position != np.floor(table.position))).any():
        raise ValueError("position must be a nonnegative integer base-pair coordinate")
    for label in ("allele_frequency", "local_ld", "imputation_quality",
                  "missingness", "temporal_coverage"):
        if ((table[label] < 0) | (table[label] > 1)).any():
            raise ValueError(f"{label} must be in [0,1]")
    if (table.recombination_rate < 0).any():
        raise ValueError("recombination_rate must be nonnegative")
    return table.set_index("locus", drop=False)


def _support(counts: pd.DataFrame) -> dict[str, frozenset[tuple[str, float]]]:
    required = {"locus", "trajectory", "time"}
    if not isinstance(counts, pd.DataFrame) or not required <= set(counts):
        raise ValueError("counts requires locus, trajectory and time")
    # The statistic validates numeric times and count data before this is called.
    return {str(locus): frozenset(zip(group.trajectory,
                                     group.time.astype(float)))
            for locus, group in counts.groupby("locus", sort=True)}


def empirical_null(
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    candidates: Sequence[str],
    *,
    family: Mapping[str, str] | None = None,
    relatedness: pd.DataFrame | None = None,
    calipers: Mapping[str, float] | None = None,
    exclusion_bp: int = 1_000_000,
    min_controls: int = 20,
    max_controls: int | None = None,
) -> EmpiricalNullResult:
    """Fit matched controls and report a finite-sample upper-tail empirical P.

    The two-sided *replication score* is ``abs(meta_slope/meta_se) *
    family_concordance``. Higher means a stronger family-discounted,
    directionally concordant change; it is not a likelihood ratio. The
    candidate's score is never used in matching. ``p_empirical`` is
    ``(1 + # control scores >= candidate score)/(m + 1)``. The +1 correction
    includes the tested locus in the exchangeable rank; ties are conservative.
    ``percentile`` is 100 times the midrank of the candidate among the m
    *controls* (strictly less plus half tied, divided by m). The P and
    percentile use different tie conventions and need not sum to 1. Small m
    limits resolution: ``p_empirical >= 1/(m+1)``.

    The entire ``candidates`` set is excluded from every other candidate's
    control pool. All counts/metadata loci are fitted, and each must have a
    defined score. Exact trajectory/time support prevents artificial changes
    in the null from missing eras or lineages; time and family labels must be
    prespecified from background information, not the tested locus.
    """
    table = _metadata(metadata)
    if (not isinstance(candidates, (list, tuple)) or not candidates
            or len(set(candidates)) != len(candidates)
            or not set(candidates) <= set(table.index)):
        raise ValueError("candidates must be a nonempty unique sequence of metadata loci")
    if not isinstance(exclusion_bp, int) or exclusion_bp < 0:
        raise ValueError("exclusion_bp must be a nonnegative integer")
    if (not isinstance(min_controls, int) or min_controls < 1
            or max_controls is not None and
            (not isinstance(max_controls, int) or max_controls < min_controls)):
        raise ValueError("require 1 <= min_controls <= max_controls (if set)")
    limits = dict(DEFAULT_CALIPERS)
    if calipers is not None:
        if set(calipers) - set(FEATURES):
            raise ValueError("unknown matching caliper")
        limits.update(calipers)
    if any(not np.isfinite(value) or value <= 0 for value in limits.values()):
        raise ValueError("calipers must be finite and positive")

    fitted = trajectory_statistic(counts, family=family, relatedness=relatedness)
    stats = fitted.summary.copy().set_index("locus", drop=False)
    if set(stats.index) != set(table.index):
        raise ValueError("metadata and fitted count loci must match exactly")
    if (stats.n_families < 2).any():
        raise ValueError("all loci require at least two prespecified trajectory families")
    stats["replication_score"] = stats.z_score.abs() * stats.family_concordance
    if not np.isfinite(stats.replication_score.to_numpy(dtype=float)).all():
        raise ValueError("all loci need a defined, finite replication score")
    support = _support(counts)
    pool = table.drop(index=list(candidates))
    result_rows = []
    match_rows = []
    for candidate in candidates:
        row = table.loc[candidate]
        eligible = pool.loc[pool.chromosome == row.chromosome].copy()
        eligible = eligible.loc[(eligible.position - row.position).abs() >= exclusion_bp]
        eligible = eligible.loc[[support[label] == support[candidate]
                                 for label in eligible.index]]
        distance = pd.Series(0.0, index=eligible.index)
        for feature in FEATURES:
            scaled = (eligible[feature] - row[feature]).abs() / limits[feature]
            eligible = eligible.loc[scaled <= 1]
            distance = distance.loc[eligible.index] + scaled.loc[eligible.index]
        if len(eligible) < min_controls:
            raise ValueError(f"{candidate}: only {len(eligible)} matched controls "
                             f"(requires {min_controls})")
        distances = pd.DataFrame({"control": eligible.index, "distance": distance.to_numpy()})
        distances = distances.sort_values(["distance", "control"], kind="stable")
        if max_controls is not None:
            distances = distances.head(max_controls)
        candidate_score = float(stats.loc[candidate, "replication_score"])
        controls = stats.loc[distances.control, "replication_score"].to_numpy(dtype=float)
        m = len(controls)
        less = int((controls < candidate_score).sum())
        equal = int((controls == candidate_score).sum())
        result_rows.append({"locus": candidate, "replication_score": candidate_score,
                            "n_controls": m, "p_empirical": (1 + m - less) / (m + 1),
                            "percentile": 100 * (less + 0.5 * equal) / m})
        for rank, (control, distance_value) in enumerate(distances.itertuples(index=False, name=None), 1):
            match_rows.append({"locus": candidate, "control": control,
                               "rank": rank, "distance": float(distance_value),
                               "control_score": float(stats.loc[control, "replication_score"])})
    return EmpiricalNullResult(pd.DataFrame(result_rows), pd.DataFrame(match_rows),
                               stats.reset_index(drop=True))
