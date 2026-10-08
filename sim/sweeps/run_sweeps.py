#!/usr/bin/env python3
"""Deterministic, time-bounded factor sweep comparing selection statistics.

For each of the six fixed-seed demographic scenarios in ``sim/scenarios`` we sweep
selection strength, uneven temporal sampling, missing-genotype rate,
trajectory-inference error and turnover rate.  Three statistics are evaluated on
the exact same sampled individuals:

* the within-trajectory replication statistic (``sim/trajectory_stat.py``) run
  the *naive* way with no relatedness matrix (working-independence sign
  concordance), labelled ``trajectory_naive``;
* the same statistic run with a supplied two-population relatedness matrix, so
  the GLS family-discount and ``n_effective`` are active, labelled
  ``trajectory_related``;
* the Akbari-style structure-corrected binomial GLMM (``sim/akbari_glmm.py``),
  labelled ``akbari``.

The sweep is one-factor-at-a-time (OFAT): a baseline cell plus, for each applicable
factor, the non-baseline levels with every other factor held at baseline.  Replicate
seeds are shared across the cells of a scenario, so factor effects are paired.

Reported quantities per (scenario, factor, level, method) are the false-positive
rate (FPR) on neutral loci, power at the designated selected locus, and effect-size
bias against the realised population frequency slope.  FPR/power are reported both
at each method's own nominal alpha and at a threshold calibrated from the neutral
``neutral_drift`` baseline (the empirical-null calibration the method documents).

Nothing here is a calibrated p-value claim: the analytic Wald p-values are working
diagnostics.  Determinism: every RNG is seeded from the run seed declared in the
scenario YAML, and the sweep is reproducible from ``(mode, runs, n_loci, seed)``.

Compute is bounded by ``runs`` and ``n_loci`` and by ``--time-budget`` seconds: the
sweep stops cleanly before starting a cell that would exceed the budget.  The tidy
``results/fpr_power.csv`` and ``results/summary.json`` are rewritten after every
completed (scenario, factor, level) cell, so an external timeout never loses work.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _path in (REPO_ROOT, REPO_ROOT / "sim", REPO_ROOT / "sim" / "scripts"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from akbari_glmm import AkbariGLMM  # noqa: E402
from simulator import DemographicSimulator, load_config  # noqa: E402
from trajectory_stat import trajectory_statistic  # noqa: E402


SCENARIO_DIR = REPO_ROOT / "sim" / "scenarios"
OUTPUT_DIR = Path(__file__).resolve().parent / "results"
FIGURE_PATH = REPO_ROOT / "fig" / "power_fpr.png"

# Fixed, declaration-order list: results must not depend on filesystem order.
SCENARIO_ORDER = (
    "neutral_drift",
    "replacement_only",
    "one_population_selection",
    "shared_selection",
    "opposing_selection",
    "selection_aligned_turnover",
)
SCENARIO_META = {
    "neutral_drift": {"has_selection": False, "has_replacement": False, "selection_pop": None},
    "replacement_only": {"has_selection": False, "has_replacement": True, "selection_pop": None},
    "one_population_selection": {"has_selection": True, "has_replacement": False, "selection_pop": "pop1"},
    "shared_selection": {"has_selection": True, "has_replacement": False, "selection_pop": "pop1"},
    "opposing_selection": {"has_selection": True, "has_replacement": False, "selection_pop": "pop1"},
    "selection_aligned_turnover": {"has_selection": True, "has_replacement": True, "selection_pop": "pop2"},
}
# Power is only meaningful where a locus is genuinely under selection.
SELECTION_SCENARIOS = tuple(
    name for name in SCENARIO_ORDER if SCENARIO_META[name]["has_selection"])
# Designated tested locus: the scenarios apply selection to locus index 0.
TEST_LOCUS_INDEX = 0

BASELINE_FACTORS = {
    "selection_strength": 1.0,
    "uneven_sampling": 1.0,
    "missing_rate": 0.0,
    "inference_error": 0.0,
    "turnover_rate": 1.0,
}
# Selection strengths are multipliers on the scenario's additive effect (0.25 per
# allele). The maximum is 1.5, not 2.0, so that opposing_selection's -0.25 effect
# stays above the simulator's -0.5 per-allele fitness-floor validity bound.
SELECTION_STRENGTH_LEVELS = (0.0, 0.5, 1.0, 1.5)
UNEVEN_SAMPLING_LEVELS = (1.0, 3.0)
MISSING_RATE_LEVELS = (0.0, 0.3)
INFERENCE_ERROR_LEVELS = (0.0, 0.15)
TURNOVER_RATE_LEVELS = (1.0, 0.5, 0.0)
ALPHA = 0.05
# Prespecified correlation of family-level slope errors between the two sampled
# populations (design choice for the sweep, not estimated at the tested locus).
RELATEDNESS_RHO = 0.5

METHODS: dict[str, dict[str, str]] = {
    "trajectory_naive": {
        "kind": "trajectory",
        "p": "traj_naive_p",
        "score": "traj_naive_score",
        "slope": "traj_naive_slope",
        "true": "true_mean_traj_slope",
    },
    "trajectory_related": {
        "kind": "trajectory",
        "p": "traj_rel_p",
        "score": "traj_rel_score",
        "slope": "traj_rel_slope",
        "true": "true_mean_traj_slope",
    },
    "akbari": {
        "kind": "akbari",
        "p": "akbari_p",
        "score": "akbari_score",
        "slope": "akbari_freq_slope",
        "true": "true_pooled_slope",
    },
}

CSV_COLUMNS = (
    "scenario", "factor", "level", "method", "n_tests",
    "fpr_nominal", "fpr_calibrated", "power_nominal", "power_calibrated",
    "effect_bias", "convergence_rate",
)


def relatedness_matrix(names: Sequence[str], rho: float = RELATEDNESS_RHO) -> pd.DataFrame:
    """Labeled positive-definite correlation matrix for the sampled populations."""
    if not 0.0 < rho < 1.0:
        raise ValueError("rho must lie strictly between 0 and 1")
    size = len(names)
    matrix = np.full((size, size), float(rho))
    np.fill_diagonal(matrix, 1.0)
    return pd.DataFrame(matrix, index=list(names), columns=list(names))


def _swept_levels(scenario: str, factor: str) -> tuple[float, ...]:
    """Applicable levels for one factor under one scenario (baseline included)."""
    meta = SCENARIO_META[scenario]
    if factor == "selection_strength":
        return SELECTION_STRENGTH_LEVELS if meta["has_selection"] else (1.0,)
    if factor == "turnover_rate":
        return TURNOVER_RATE_LEVELS if meta["has_replacement"] else (1.0,)
    if factor == "uneven_sampling":
        return UNEVEN_SAMPLING_LEVELS
    if factor == "missing_rate":
        return MISSING_RATE_LEVELS
    if factor == "inference_error":
        return INFERENCE_ERROR_LEVELS
    raise KeyError(factor)


def build_cells(scenario: str) -> list[dict[str, Any]]:
    """One-factor-at-a-time cells: baseline plus each non-baseline factor level."""
    cells = [{"factor": "baseline", "level": 1.0, "factors": dict(BASELINE_FACTORS)}]
    for factor in BASELINE_FACTORS:
        for level in _swept_levels(scenario, factor):
            if level == BASELINE_FACTORS[factor]:
                continue
            factors = dict(BASELINE_FACTORS)
            factors[factor] = float(level)
            cells.append({"factor": factor, "level": float(level), "factors": factors})
    return cells


def perturb_config(config: dict[str, Any], factors: dict[str, float]) -> dict[str, Any]:
    """Return a deep copy of the scenario config with the sampled factors applied.

    ``selection_strength`` scales every additive selection effect; ``turnover_rate``
    scales every replacement pulse fraction.  Uneven sampling, missing genotypes and
    trajectory-inference error are observation-level perturbations applied later, so
    they do not appear here.  A scenario without a feature treats the factor as a
    no-op (recorded as not applicable by the cell grid).
    """
    updated = copy.deepcopy(config)
    multiplier = float(factors["selection_strength"])
    selection = updated.get("selection") or {}
    updated["selection"] = {
        pop: {locus: float(effect) * multiplier for locus, effect in loci.items()}
        for pop, loci in selection.items()
    }
    turnover = float(factors["turnover_rate"])
    updated["replacement"] = [
        dict(event, fraction=float(event["fraction"]) * turnover)
        for event in updated.get("replacement", [])
    ]
    # Loci are treated as unlinked so the background GRM is a pure structure term,
    # not a local-LD echo of the tested locus.  This is a sweep design choice.
    updated["recombination_rate"] = 0.5
    return updated


class RecordingSimulator(DemographicSimulator):
    """Simulator that also records the realised population allele frequencies.

    Used only to obtain the true per-generation frequencies at sampled generations,
    enabling an honest effect-size bias against the realised trajectory rather than
    against a design coefficient.
    """

    def __init__(self, config: dict[str, Any]):
        super().__init__(config)
        self.frequencies: list[tuple[int, str, int, float]] = []

    def _sample(self, generation: int, rows: list[dict[str, Any]]) -> None:
        super()._sample(generation, rows)
        for name in self.names:
            pool = self.populations[name]
            counts = pool.sum(axis=(1, 2))
            totals = 2 * pool.shape[0]
            for locus, count in enumerate(counts):
                self.frequencies.append((generation, name, int(locus), float(count / totals)))


def _ols_slope(x: np.ndarray, y: np.ndarray) -> float:
    """Simple least-squares slope of ``y`` on ``x`` (equal weights)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xc = x - x.mean()
    denom = float(xc @ xc)
    if denom <= 0:
        return float("nan")
    return float(xc @ y / denom)


def _true_slopes(sim: RecordingSimulator) -> tuple[pd.DataFrame, dict[tuple[int, str], float], pd.DataFrame]:
    """Per-trajectory and pooled realised frequency slopes per locus (per year BP)."""
    frame = pd.DataFrame(sim.frequencies, columns=["generation", "trajectory", "locus", "frequency"])
    years = (sim.generations - frame.generation) * sim.config.get("generation_years", 25)
    frame["years_bp"] = years
    per_traj = {}
    for (locus, trajectory), group in frame.groupby(["locus", "trajectory"], sort=True):
        per_traj[(locus, trajectory)] = _ols_slope(group.years_bp.to_numpy(), group.frequency.to_numpy())
    pooled = []
    for locus, group in frame.groupby("locus", sort=True):
        by_gen = group.groupby("generation", as_index=False).frequency.mean()
        years_gen = (sim.generations - by_gen.generation) * sim.config.get("generation_years", 25)
        pooled.append(_ols_slope(years_gen.to_numpy(), by_gen.frequency.to_numpy()))
    table = pd.DataFrame(
        [{"locus": f"locus_{l}",
          "true_pooled_slope": pooled[l],
          "true_mean_traj_slope": float(np.mean([per_traj[(l, t)] for t in sim.names]))}
         for l in range(sim.n_loci)]
    )
    return table, per_traj, frame


def _apply_observation_perturbations(
    rows: list[dict[str, Any]],
    sim: RecordingSimulator,
    factors: dict[str, float],
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (genotypes, inferred trajectories, times, keep_mask, missing_mask).

    Uneven sampling thins internal generations by ``1/uneven_sampling``; the endpoint
    generations are retained to keep the temporal span.  Missing genotypes are drawn
    per individual and locus.  Trajectory-inference error relabels a fraction of
    individuals with a different (here: the other) inferred trajectory.
    """
    rows = list(rows)
    genotypes = np.array([r["genotype"] for r in rows], dtype=float)
    times = np.array([r["years_bp"] for r in rows], dtype=float)
    generations = np.array([r["generation"] for r in rows], dtype=int)
    trajectories = np.array([r["trajectory"] for r in rows], dtype=object)

    keep = np.ones(len(rows), dtype=bool)
    uneven = float(factors["uneven_sampling"])
    if uneven > 1.0:
        internal = (generations != generations.min()) & (generations != generations.max())
        keep &= (~internal) | (rng.random(len(rows)) < 1.0 / uneven)

    inference_error = float(factors["inference_error"])
    if inference_error > 0:
        names = list(sim.names)
        order = {name: index for index, name in enumerate(names)}
        swap = rng.random(len(rows)) < inference_error
        next_index = np.array([(order[name] + 1) % len(names) for name in trajectories])
        trajectories = np.where(swap, np.array(names, dtype=object)[next_index], trajectories)

    missing_rate = float(factors["missing_rate"])
    missing = rng.random(genotypes.shape) < missing_rate
    return genotypes, trajectories, times, keep, missing


def _trajectory_counts(
    genotypes: np.ndarray, trajectories: np.ndarray, times: np.ndarray,
    keep: np.ndarray, missing: np.ndarray, locus_labels: np.ndarray,
) -> pd.DataFrame:
    """Build the long locus/trajectory/time counts fed to ``trajectory_statistic``."""
    observed_alt = genotypes * (~missing)
    observed_total = 2.0 * (~missing)
    rows_kept = np.flatnonzero(keep)
    long = pd.DataFrame({
        "locus": np.tile(locus_labels, len(rows_kept)),
        "trajectory": np.repeat(trajectories[rows_kept], len(locus_labels)),
        "time": np.repeat(times[rows_kept], len(locus_labels)),
        "alt_count": observed_alt[rows_kept].ravel(),
        "total_count": observed_total[rows_kept].ravel(),
    })
    counts = (long.groupby(["locus", "trajectory", "time"], as_index=False)[["alt_count", "total_count"]]
              .sum())
    counts = counts[counts.total_count > 0].reset_index(drop=True)
    counts[["alt_count", "total_count"]] = counts[["alt_count", "total_count"]].astype(int)
    # Drop any (locus, trajectory) with fewer than two observed times; the statistic
    # requires at least two distinct times and otherwise fails closed.
    informative = counts.groupby(["locus", "trajectory"]).time.transform("nunique") >= 2
    return counts[informative].reset_index(drop=True)


def run_replicate(
    config: dict[str, Any],
    factors: dict[str, float],
    seed: int,
    akbari_loci: int,
    relatedness: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Simulate one replicate and return one tidy row per tested locus."""
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), 0xC0FFEE]))
    sim = RecordingSimulator(config)
    sim.seed = int(seed)
    sim.rng = np.random.default_rng(int(seed))
    result = sim.run()
    rows = result["samples"]

    genotypes, trajectories, times, keep, missing = _apply_observation_perturbations(
        rows, sim, factors, rng)
    n_loci = sim.n_loci
    locus_labels = np.array([f"locus_{l}" for l in range(n_loci)])

    # ---- within-trajectory replication statistic, two relatedness treatments ----
    counts = _trajectory_counts(genotypes, trajectories, times, keep, missing, locus_labels)
    naive = trajectory_statistic(counts).summary.set_index("locus")
    related = trajectory_statistic(counts, relatedness=relatedness).summary.set_index("locus")

    # ---- Akbari-style GLMM -----------------------------------------------------
    akbari_genotypes = np.where(missing, np.nan, genotypes)
    akbari_rows: dict[str, dict[str, float]] = {}
    n_akbari = min(akbari_loci, n_loci)
    for l in range(n_akbari):
        label = f"locus_{l}"
        observed = keep & ~missing[:, l]
        if observed.sum() < 6:
            akbari_rows[label] = {"akbari_z": np.nan, "akbari_p": np.nan,
                                  "akbari_freq_slope": np.nan, "akbari_converged": 0.0}
            continue
        background = np.ones(n_loci, dtype=bool)
        background[l] = False
        k = AkbariGLMM.compute_grm(akbari_genotypes[observed][:, background])
        data = pd.DataFrame({
            "time_bp": times[observed],
            "alt_count": genotypes[observed, l],
            "total_count": np.full(int(observed.sum()), 2.0),
        })
        try:
            fit = AkbariGLMM().fit(data, K=k)
        except (RuntimeError, ValueError):
            akbari_rows[label] = {"akbari_z": np.nan, "akbari_p": np.nan,
                                  "akbari_freq_slope": np.nan, "akbari_converged": 0.0}
            continue
        pbar = float(genotypes[observed, l].sum() / data.total_count.sum())
        akbari_rows[label] = {
            "akbari_z": fit["z_score"], "akbari_p": fit["p_value"],
            "akbari_freq_slope": pbar * (1 - pbar) * fit["beta_time"],
            "akbari_converged": 1.0,
        }

    truth, _per_traj, _frame = _true_slopes(sim)
    truth = truth.set_index("locus")

    records = []
    for locus in sorted(naive.index):
        index = int(locus.split("_")[1])
        naive_row = naive.loc[locus]
        related_row = related.loc[locus] if locus in related.index else naive_row
        a = akbari_rows.get(locus, {"akbari_z": np.nan, "akbari_p": np.nan,
                                    "akbari_freq_slope": np.nan, "akbari_converged": 0.0})

        def _score(row: pd.Series) -> float:
            concordance = float(row.family_concordance)
            if not np.isfinite(concordance):
                return float("nan")
            return abs(float(row.z_score)) * concordance

        records.append({
            "locus": locus,
            "is_selected": index == TEST_LOCUS_INDEX,
            "traj_naive_z": float(naive_row.z_score),
            "traj_naive_p": float(naive_row.p_value),
            "traj_naive_score": _score(naive_row),
            "traj_naive_slope": float(naive_row.meta_slope),
            "traj_rel_z": float(related_row.z_score),
            "traj_rel_p": float(related_row.p_value),
            "traj_rel_score": _score(related_row),
            "traj_rel_slope": float(related_row.meta_slope),
            "traj_n_effective": float(related_row.n_effective),
            "true_pooled_slope": float(truth.loc[locus, "true_pooled_slope"]),
            "true_mean_traj_slope": float(truth.loc[locus, "true_mean_traj_slope"]),
            "akbari_z": a["akbari_z"], "akbari_p": a["akbari_p"],
            "akbari_score": abs(a["akbari_z"]) if np.isfinite(a["akbari_z"]) else np.nan,
            "akbari_freq_slope": a["akbari_freq_slope"],
            "akbari_converged": a["akbari_converged"],
        })
    frame = pd.DataFrame.from_records(records)
    frame["replicate"] = int(seed)
    return frame


def _threshold(reference: pd.DataFrame, column: str) -> float:
    values = reference[column].to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    return float(np.quantile(values, 1.0 - ALPHA))


def calibration_thresholds(locus_rows: pd.DataFrame) -> dict[str, float]:
    """Per-method 95th-percentile neutral scores from the ``neutral_drift`` baseline.

    The baseline cell is therefore in-sample for that calibration; every other cell
    is an out-of-sample transfer test.
    """
    reference = locus_rows[(locus_rows.scenario == "neutral_drift")
                           & (locus_rows.factor == "baseline")
                           & (~locus_rows.is_selected)]
    thresholds = {}
    for method, spec in METHODS.items():
        rows = reference
        if spec["kind"] == "akbari":
            rows = rows[rows.akbari_converged == 1.0]
        thresholds[method] = _threshold(rows, spec["score"])
    return thresholds


def _cell_metric_rows(cell: pd.DataFrame, scenario: str, factor: str, level: float,
                      thresholds: dict[str, float]) -> list[dict[str, Any]]:
    neutral = cell[~cell.is_selected]
    # Power/bias are only defined where a locus is genuinely under selection.
    has_selection = SCENARIO_META[scenario]["has_selection"]
    selected = cell[cell.is_selected] if has_selection else cell.iloc[0:0]
    rows = []
    for method, spec in METHODS.items():
        base = neutral
        conv_cell = cell
        if spec["kind"] == "akbari":
            base = neutral[neutral.akbari_converged == 1.0]
            conv_cell = cell[cell.akbari_converged == 1.0]
        # --- calibration / convergence rates ---
        if spec["kind"] == "akbari":
            convergence = float(cell.akbari_converged.mean()) if len(cell) else np.nan
        else:
            score_all = cell[spec["score"]].to_numpy(dtype=float)
            convergence = float(np.isfinite(score_all).mean()) if len(cell) else np.nan
        fpr_nominal = fpr_calibrated = power_nominal = power_calibrated = np.nan
        effect_bias = np.nan
        # --- FPR on neutral loci ---
        valid_p = base[np.isfinite(base[spec["p"]].to_numpy(dtype=float))]
        if len(valid_p):
            fpr_nominal = float((valid_p[spec["p"]] < ALPHA).mean())
        score_valid = base[np.isfinite(base[spec["score"]].to_numpy(dtype=float))]
        threshold = thresholds.get(method, float("nan"))
        n_tests = int(len(score_valid))
        if len(score_valid) and np.isfinite(threshold):
            fpr_calibrated = float((score_valid[spec["score"]] > threshold).mean())
        # --- power and bias at the selected locus (selection scenarios only) ---
        if len(selected):
            sel = selected[selected.akbari_converged == 1.0] if spec["kind"] == "akbari" else selected
            sel_p = sel[np.isfinite(sel[spec["p"]].to_numpy(dtype=float))]
            if len(sel_p):
                power_nominal = float((sel_p[spec["p"]] < ALPHA).mean())
            sel_score = sel[np.isfinite(sel[spec["score"]].to_numpy(dtype=float))]
            if len(sel_score) and np.isfinite(threshold):
                power_calibrated = float((sel_score[spec["score"]] > threshold).mean())
            bias = (selected[spec["slope"]] - selected[spec["true"]]).to_numpy(dtype=float)
            bias = bias[np.isfinite(bias)]
            if bias.size:
                effect_bias = float(bias.mean())
        rows.append({
            "scenario": scenario, "factor": factor, "level": float(level),
            "method": method, "n_tests": n_tests,
            "fpr_nominal": fpr_nominal, "fpr_calibrated": fpr_calibrated,
            "power_nominal": power_nominal, "power_calibrated": power_calibrated,
            "effect_bias": effect_bias, "convergence_rate": convergence,
        })
    return rows


def compute_metrics(locus_rows: pd.DataFrame,
                    thresholds: dict[str, float] | None = None) -> tuple[pd.DataFrame, dict[str, float]]:
    """Aggregate FPR, power and effect-size bias per (scenario, factor, level, method)."""
    if thresholds is None:
        thresholds = calibration_thresholds(locus_rows)
    rows: list[dict[str, Any]] = []
    for (scenario, factor, level), cell in locus_rows.groupby(
            ["scenario", "factor", "level"], sort=True):
        rows.extend(_cell_metric_rows(cell, scenario, factor, float(level), thresholds))
    return pd.DataFrame(rows, columns=list(CSV_COLUMNS)), thresholds


def cell_grid(scenarios: Sequence[str]) -> list[tuple[str, dict[str, Any]]]:
    """Ordered (scenario, cell) pairs; ``neutral_drift`` baseline leads the list."""
    grid = []
    for scenario in scenarios:
        if scenario not in SCENARIO_META:
            raise ValueError(f"unknown scenario {scenario!r}")
        for cell in build_cells(scenario):
            grid.append((scenario, cell))
    return grid


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def _atomic_write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False)
    os.replace(tmp, path)


def _claims(metrics: pd.DataFrame) -> list[str]:
    """Short, deterministic claims read straight off the aggregated CSV."""
    claims: list[str] = []

    def _row(scenario: str, factor: str, method: str) -> pd.Series | None:
        hit = metrics[(metrics.scenario == scenario) & (metrics.factor == factor)
                      & (metrics.method == method)]
        return None if hit.empty else hit.iloc[0]

    rep_naive = _row("replacement_only", "baseline", "trajectory_naive")
    rep_rel = _row("replacement_only", "baseline", "trajectory_related")
    neutral_naive = _row("neutral_drift", "baseline", "trajectory_naive")
    if rep_naive is not None and rep_rel is not None and neutral_naive is not None:
        claims.append(
            "observation: replacement_only calibrated FPR is "
            f"{rep_naive.fpr_calibrated:.3f} for trajectory_naive vs "
            f"{neutral_naive.fpr_calibrated:.3f} under neutral_drift (nominal alpha "
            f"{ALPHA}); supplying the relatedness matrix changes it to "
            f"{rep_rel.fpr_calibrated:.3f}, a repair of "
            f"{rep_naive.fpr_calibrated - rep_rel.fpr_calibrated:.3f}.")
    for scenario in SELECTION_SCENARIOS:
        selection_rows = metrics[(metrics.scenario == scenario) & (metrics.factor == "baseline")]
        if selection_rows.empty:
            continue
        best = selection_rows.sort_values("power_calibrated", ascending=False)
        top = best.iloc[0]
        akbari = selection_rows[selection_rows.method == "akbari"]
        if not akbari.empty:
            gap = float(top.power_calibrated) - float(akbari.iloc[0].power_calibrated)
            claims.append(
                f"inference: {scenario} baseline calibrated power is highest for "
                f"{top.method} ({top.power_calibrated:.3f}); akbari is "
                f"{akbari.iloc[0].power_calibrated:.3f} (gap {gap:+.3f}).")
    return claims


def _summary_payload(design: dict[str, Any], metrics: pd.DataFrame) -> dict[str, Any]:
    payload = dict(design)
    payload["claims"] = _claims(metrics)
    payload["n_cells_scored"] = int(len(metrics.drop_duplicates(["scenario", "factor", "level"])))
    return payload


def run_sweep(
    scenarios: Sequence[str],
    *,
    runs: int,
    n_loci: int,
    samples_per_population: int,
    label: str,
    akbari_loci: int = 6,
    rho: float = RELATEDNESS_RHO,
    output_dir: Path = OUTPUT_DIR,
    figure_path: Path | None = FIGURE_PATH,
    time_budget: float | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run the OFAT sweep, writing CSV+JSON incrementally after every cell.

    Returns the tidy metrics frame and the design/summary payload.  If
    ``time_budget`` seconds elapse, the loop stops before the next cell and writes
    whatever has completed; ``design['partial']`` then records the truncation.
    """
    if not runs >= 1 or n_loci < 4 or samples_per_population < 2:
        raise ValueError("require runs >= 1, n_loci >= 4, samples_per_population >= 2")
    if akbari_loci < 1:
        raise ValueError("akbari_loci must be at least 1")
    grid = cell_grid(scenarios)
    scenario_index = {name: i for i, name in enumerate(SCENARIO_ORDER)}
    started = time.time()
    thresholds: dict[str, float] | None = None
    metric_rows: list[dict[str, Any]] = []
    completed: list[str] = []
    csv_path = output_dir / "fpr_power.csv"
    summary_path = output_dir / "summary.json"

    design: dict[str, Any] = {
        "label": label,
        "runs": int(runs),
        "n_loci": int(n_loci),
        "samples_per_population": int(samples_per_population),
        "akbari_loci": int(min(akbari_loci, n_loci)),
        "alpha": ALPHA,
        "relatedness_rho": float(rho),
        "scenarios": list(scenarios),
        "baseline_factors": BASELINE_FACTORS,
        "selection_strength_levels": list(SELECTION_STRENGTH_LEVELS),
        "uneven_sampling_levels": list(UNEVEN_SAMPLING_LEVELS),
        "missing_rate_levels": list(MISSING_RATE_LEVELS),
        "inference_error_levels": list(INFERENCE_ERROR_LEVELS),
        "turnover_rate_levels": list(TURNOVER_RATE_LEVELS),
        "scenario_seeds": {name: int(load_config(SCENARIO_DIR / f"{name}.yaml")["seed"])
                           for name in scenarios},
        "replicate_seed_formula": "scenario_seed * 10000 + replicate",
        "perturbation_seed_formula": "SeedSequence([replicate_seed, 0xC0FFEE])",
        "n_cells_planned": len(grid),
        "calibration": "neutral_drift baseline neutral loci, 95th percentile score",
    }
    design["wall_seconds"] = 0.0
    design["partial"] = False

    for scenario, cell in grid:
        if time_budget is not None and (time.time() - started) > time_budget:
            design["partial"] = True
            design["stop_reason"] = "time budget reached before starting next cell"
            break
        config = load_config(SCENARIO_DIR / f"{scenario}.yaml")
        config["n_loci"] = int(n_loci)
        config["samples_per_population"] = int(samples_per_population)
        config = perturb_config(config, cell["factors"])
        relatedness = relatedness_matrix(list(config["populations"]), rho)
        frames = []
        for replicate in range(runs):
            seed = int(config["seed"]) * 10_000 + replicate
            frame = run_replicate(config, cell["factors"], seed,
                                  akbari_loci=akbari_loci, relatedness=relatedness)
            frame["scenario"] = scenario
            frame["factor"] = cell["factor"]
            frame["level"] = float(cell["level"])
            frame["scenario_index"] = scenario_index[scenario]
            frames.append(frame)
        cell_rows = pd.concat(frames, ignore_index=True)
        if thresholds is None and scenario == "neutral_drift" and cell["factor"] == "baseline":
            thresholds = calibration_thresholds(cell_rows)
            design["calibration_thresholds"] = thresholds
        if thresholds is None:
            # No neutral reference yet (e.g. neutral_drift excluded): score nominal only.
            thresholds = {method: float("nan") for method in METHODS}
        metric_rows.extend(_cell_metric_rows(cell_rows, scenario, cell["factor"],
                                             float(cell["level"]), thresholds))
        completed.append(f"{scenario}:{cell['factor']}={cell['level']:g}")
        metrics = pd.DataFrame(metric_rows, columns=list(CSV_COLUMNS))
        design["wall_seconds"] = round(time.time() - started, 2)
        design["cells_completed"] = list(completed)
        _atomic_write_csv(metrics, csv_path)
        _atomic_write_text(summary_path,
                           json.dumps(_summary_payload(design, metrics), indent=2, sort_keys=True) + "\n")

    metrics = pd.DataFrame(metric_rows, columns=list(CSV_COLUMNS))
    design["wall_seconds"] = round(time.time() - started, 2)
    design["cells_completed"] = list(completed)
    _atomic_write_csv(metrics, csv_path)
    _atomic_write_text(summary_path,
                       json.dumps(_summary_payload(design, metrics), indent=2, sort_keys=True) + "\n")
    if figure_path is not None and not metrics.empty:
        make_figure(csv_path, figure_path)
    return metrics, design


def make_figure(csv_path: Path, figure_path: Path) -> None:
    """Regenerate the FPR/power heatmaps from the tidy ``fpr_power.csv``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = pd.read_csv(csv_path)
    metrics["cell"] = metrics.apply(
        lambda r: r.scenario if r.factor == "baseline" else f"{r.scenario}\n{r.factor}={r.level:g}",
        axis=1)
    method_order = list(METHODS)
    method_short = {"trajectory_naive": "traj naive", "trajectory_related": "traj related",
                    "akbari": "akbari"}

    fpr = (metrics.pivot_table(index="cell", columns="method", values="fpr_calibrated",
                               aggfunc="first").reindex(columns=method_order))
    power_cells = metrics[metrics.scenario.isin(SELECTION_SCENARIOS)]
    power = (power_cells.pivot_table(index="cell", columns="method", values="power_calibrated",
                                     aggfunc="first").reindex(columns=method_order))

    fig, axes = plt.subplots(1, 2, figsize=(13, max(6, 0.34 * len(fpr))),
                             gridspec_kw={"width_ratios": [1.15, 1.0]})
    for ax, data, title, vmax, cmap in (
        (axes[0], fpr, f"Calibrated FPR (neutral loci)\nnominal alpha={ALPHA}", 0.5, "magma_r"),
        (axes[1], power, "Calibrated power (selected locus)", 1.0, "viridis"),
    ):
        matrix = data.to_numpy(dtype=float)
        image = ax.imshow(matrix, aspect="auto", vmin=0.0, vmax=vmax, cmap=cmap)
        ax.set_xticks(range(len(method_order)))
        ax.set_xticklabels([method_short[m] for m in method_order], rotation=30, ha="right")
        ax.set_yticks(range(len(data.index)))
        ax.set_yticklabels(data.index, fontsize=6)
        ax.set_title(title)
        for row in range(matrix.shape[0]):
            for col in range(matrix.shape[1]):
                value = matrix[row, col]
                if np.isfinite(value):
                    ax.text(col, row, f"{value:.2f}", ha="center", va="center", fontsize=6,
                            color="white" if value > 0.6 * vmax else "black")
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Replication statistic (naive / relatedness-discounted) vs Akbari GLMM")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=150)
    plt.close(fig)


QUICK_DEFAULTS = {"runs": 6, "n_loci": 12, "samples_per_population": 8, "akbari_loci": 3}
FULL_DEFAULTS = {"runs": 40, "n_loci": 24, "samples_per_population": 10, "akbari_loci": 6}


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--quick", action="store_true",
                      help="small bounded run used for validation (<90s)")
    mode.add_argument("--full", action="store_true",
                      help="report run with the full replicate/locus budget")
    parser.add_argument("--runs", type=int, default=None, help="replicates per cell")
    parser.add_argument("--n-loci", type=int, default=None, help="simulated loci per replicate")
    parser.add_argument("--samples-per-population", type=int, default=None)
    parser.add_argument("--akbari-loci", type=int, default=None,
                        help="loci (including the tested locus) fitted by the GLMM")
    parser.add_argument("--rho", type=float, default=RELATEDNESS_RHO,
                        help="prespecified family slope-error correlation")
    parser.add_argument("--time-budget", type=float, default=None,
                        help="wall-clock seconds before stopping between cells")
    parser.add_argument("--scenarios", nargs="+", default=list(SCENARIO_ORDER))
    parser.add_argument("--output", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--figure", type=Path, default=FIGURE_PATH)
    parser.add_argument("--no-figure", action="store_true")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    quick = bool(args.quick)
    defaults = QUICK_DEFAULTS if quick else FULL_DEFAULTS
    label = "quick" if quick else "full"
    runs = args.runs if args.runs is not None else defaults["runs"]
    n_loci = args.n_loci if args.n_loci is not None else defaults["n_loci"]
    samples = (args.samples_per_population if args.samples_per_population is not None
               else defaults["samples_per_population"])
    akbari_loci = args.akbari_loci if args.akbari_loci is not None else defaults["akbari_loci"]
    budget = args.time_budget
    if budget is None:
        budget = 75.0 if quick else 660.0
    figure_path = None if args.no_figure else args.figure
    metrics, design = run_sweep(
        args.scenarios, runs=runs, n_loci=n_loci, samples_per_population=samples,
        label=label, akbari_loci=akbari_loci, rho=args.rho, output_dir=args.output,
        figure_path=figure_path, time_budget=budget)
    print(f"[{label}] {len(design['cells_completed'])}/{design['n_cells_planned']} cells x "
          f"{runs} replicates in {design['wall_seconds']}s; rows={len(metrics)}; "
          f"partial={design['partial']}")
    print(json.dumps(design.get("calibration_thresholds", {}), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
