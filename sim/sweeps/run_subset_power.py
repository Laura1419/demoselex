#!/usr/bin/env python3
"""Subset-selection power/FPR sweep: max-over-subsets vs the structure-corrected GLMM.

Design
------
``sim/scenarios/subset_selection.yaml`` simulates K = 6 independent trajectories
(populations).  For each cell of the grid (f, strength) the tested locus (index 0)
is under positive selection in the first ``n_selected = round(f*K)`` trajectories
and drifts neutrally in the rest; ``strength`` scales the baseline additive effect
0.25/allele.  The tested locus is index 0 and every other locus is a neutral
control locus.

Three statistics are scored on the same sampled individuals:

* ``akbari`` - the Akbari-style structure-corrected binomial time/GRM GLMM
  (``sim/akbari_glmm.py``), the "structure-corrected GLMM" reference;
* ``traj_all`` - the locus meta statistic of ``sim/trajectory_stat.py`` over ALL
  trajectories, with working independence (identity relatedness), i.e. the plain
  concordance/meta statistic.  ``score = |z_score| * family_concordance``;
* ``traj_max_subset`` - the max-over-subsets device.

max-over-subsets definition
---------------------------
``analysis/subset_regime.py`` (and ``docs/subset-regime-prereg.md``) do NOT exist
in this tree, so the device is defined here, from the same per-trajectory slopes
that ``trajectory_stat.py`` produces, and the definition is stated explicitly:

    For a locus, let (beta_i, se_i) be the weighted frequency slope and binomial
    SE of trajectory i.  For every subset S of the K trajectories with |S| >= 2,
    form the working-independence inverse-variance meta slope on S,
        meta_S = sum_{i in S} beta_i/se_i^2 / sum_{i in S} 1/se_i^2,
        z_S    = meta_S * sqrt(sum_{i in S} 1/se_i^2),
    and its sign concordance conc_S = max(n+, n-)/(n+ + n-) over S (undefined,
    hence skipped, when every slope is exactly 0).  The device statistic is
        score = max_{S: |S|>=2} |z_S| * conc_S,
    and ``max_subset_size`` records |argmax S|.

This is exactly ``trajectory_statistic``'s meta/GLS formula run on a subset with an
identity relatedness matrix (no family collapse; each trajectory is its own family),
so the individual subset scores are numerically the same quantity the all-trajectory
statistic uses.  Singletons are excluded because a one-slope concordance is
undefined, so the device can only fire when at least two trajectories agree.

Calibration
-----------
Matched-null calibration follows ``sim/sweeps/run_sweeps.py``: the 95th percentile
of each statistic's score on the neutral loci of the neutral regime is used as the
threshold.  The neutral reference is the same K = 6 population set with selection
removed, so it is matched in K and sampling to the tested cells.  The named
``neutral_drift`` and ``replacement_only`` scenarios have only two trajectories, so
they cannot exercise subset search; their demographic regimes are therefore lifted
onto the same K = 6 set to measure the device's false-positive cost:

* ``neutral_drift``  - K = 6, no selection, no replacement (independent drift);
* ``replacement_only`` - K = 6, no selection, plus replacement_only's pulse
  (pop1 <- pop2 at generation 12, fraction 0.85).

``fpr_naive_device`` is the device's FPR when its maximum is judged against the
single all-trajectory 95th-percentile threshold (no multiplicity correction), the
cost you pay for scanning subsets without recalibrating the maximum.

Definitions and negative results are first-class; nothing here is a calibrated
genomic p-value claim.  Output is written incrementally after every cell so an
external timeout never loses work.  Determinism: replicate seed = scenario_seed *
10000 + replicate.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _path in (REPO_ROOT, REPO_ROOT / "sim", REPO_ROOT / "sim" / "scripts",
              REPO_ROOT / "sim" / "sweeps"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from akbari_glmm import AkbariGLMM  # noqa: E402
from run_sweeps import (  # noqa: E402
    ALPHA, BASELINE_FACTORS, RecordingSimulator, _apply_observation_perturbations,
    _atomic_write_csv, _atomic_write_text, _trajectory_counts,
)
from simulator import load_config  # noqa: E402
from trajectory_stat import trajectory_statistic  # noqa: E402


SCENARIO_PATH = REPO_ROOT / "sim" / "scenarios" / "subset_selection.yaml"
OUTPUT_DIR = Path(__file__).resolve().parent / "results"
CSV_PATH = OUTPUT_DIR / "subset_power.csv"
SUMMARY_PATH = OUTPUT_DIR / "summary.json"
FIGURE_PATH = REPO_ROOT / "fig" / "subset_power.png"

TEST_LOCUS_INDEX = 0
MIN_SUBSET_SIZE = 2
# Baseline additive effect in subset_selection.yaml; ``strength`` is a multiplier.
BASELINE_EFFECT = 0.25
F_LEVELS = (0.2, 0.4, 0.6, 1.0)
STRENGTH_LEVELS = (0.5, 1.0, 2.0)
METHODS = ("akbari", "traj_all", "traj_max_subset")
# Demographic regime lifted onto the K = 6 population set (see module docstring).
REPLACEMENT_EVENT = {"generation": 12, "recipient": "pop1", "source": "pop2", "fraction": 0.85}

CSV_COLUMNS = (
    "scenario", "f", "strength", "n_trajectories", "n_selected", "method", "n_tests",
    "fpr_nominal", "fpr_calibrated", "fpr_naive_device",
    "power_nominal", "power_calibrated", "power_naive_device",
    "mean_best_subset_size", "convergence_rate",
)


# --------------------------------------------------------------------------- design
def subset_config(base: dict[str, Any], n_loci: int, samples: int,
                  n_selected: int, strength: float,
                  replace: bool = False) -> dict[str, Any]:
    """Return the subset_selection config with ``n_selected`` trajectories selected.

    The first ``n_selected`` populations in declaration order are selected at
    ``BASELINE_EFFECT * strength`` on locus 0; the rest are left to drift.  Loci are
    treated as unlinked (recombination 0.5) so the background GRM is a pure structure
    term, matching the ``run_sweeps`` design choice.
    """
    updated = copy.deepcopy(base)
    updated["n_loci"] = int(n_loci)
    updated["samples_per_population"] = int(samples)
    updated["recombination_rate"] = 0.5
    names = list(updated["populations"])
    if not 0 <= n_selected <= len(names):
        raise ValueError("n_selected must lie within the trajectory set")
    updated["selection"] = {name: {"0": BASELINE_EFFECT * float(strength)}
                            for name in names[:n_selected]}
    updated["replacement"] = [dict(REPLACEMENT_EVENT)] if replace else []
    return updated


def _fraction_to_count(fraction: float, k: int) -> int:
    """Round f*K to the nearest trajectory count (deterministic, half-up)."""
    return int(np.floor(float(fraction) * k + 0.5))


def _subset_combinations(size: int) -> list[tuple[int, ...]]:
    return [combo for r in range(MIN_SUBSET_SIZE, size + 1)
            for combo in itertools.combinations(range(size), r)]


def _max_subset_score(beta: np.ndarray, se: np.ndarray) -> tuple[float, int]:
    """Max over subsets of |meta z| * sign concordance (working independence)."""
    weights = 1.0 / np.square(se)
    best = float("nan")
    best_size = 0
    for combo in _subset_combinations(len(beta)):
        idx = list(combo)
        b = beta[idx]
        w = weights[idx]
        info = float(w.sum())
        z = float(b @ w) / np.sqrt(info)
        positive = int(np.count_nonzero(b > 0))
        negative = int(np.count_nonzero(b < 0))
        total = positive + negative
        if total == 0:
            continue
        score = abs(z) * (max(positive, negative) / total)
        if not np.isfinite(best) or score > best:
            best = score
            best_size = len(combo)
    return best, best_size


# ------------------------------------------------------------------------ one replicate
def run_replicate(config: dict[str, Any], seed: int, akbari_loci: int) -> pd.DataFrame:
    """Simulate one replicate and return one tidy row per locus with all statistics."""
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), 0xC0FFEE]))
    sim = RecordingSimulator(config)
    sim.seed = int(seed)
    sim.rng = np.random.default_rng(int(seed))
    result = sim.run()
    rows = result["samples"]

    genotypes, trajectories, times, keep, missing = _apply_observation_perturbations(
        rows, sim, BASELINE_FACTORS, rng)
    n_loci = sim.n_loci
    locus_labels = np.array([f"locus_{l}" for l in range(n_loci)])

    counts = _trajectory_counts(genotypes, trajectories, times, keep, missing, locus_labels)
    traj = trajectory_statistic(counts)
    summary = traj.summary.set_index("locus")
    slopes = traj.slopes

    subset_stats: dict[str, tuple[float, int]] = {}
    for locus, group in slopes.groupby("locus", sort=True):
        group = group.sort_values("trajectory")
        subset_stats[locus] = _max_subset_score(
            group.slope.to_numpy(dtype=float), group.se.to_numpy(dtype=float))

    # ---- Akbari-style GLMM (same wiring as run_sweeps.run_replicate) -------------
    akbari_genotypes = np.where(missing, np.nan, genotypes)
    akbari_rows: dict[str, dict[str, float]] = {}
    for l in range(min(akbari_loci, n_loci)):
        label = f"locus_{l}"
        observed = keep & ~missing[:, l]
        if observed.sum() < 6:
            akbari_rows[label] = {"akbari_z": np.nan, "akbari_p": np.nan,
                                  "akbari_converged": 0.0}
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
                                  "akbari_converged": 0.0}
            continue
        akbari_rows[label] = {"akbari_z": fit["z_score"], "akbari_p": fit["p_value"],
                              "akbari_converged": 1.0}

    records = []
    for locus in sorted(summary.index):
        index = int(locus.split("_")[1])
        row = summary.loc[locus]
        concordance = float(row.family_concordance)
        traj_score = (abs(float(row.z_score)) * concordance
                      if np.isfinite(concordance) else float("nan"))
        best, best_size = subset_stats.get(locus, (float("nan"), 0))
        akbari = akbari_rows.get(locus)
        if akbari is None:
            # Locus never fitted by the GLMM (outside the akbari_loci budget):
            # leave its columns NaN so convergence is averaged over fitted loci only.
            akbari = {"akbari_z": np.nan, "akbari_p": np.nan, "akbari_converged": np.nan}
        records.append({
            "locus": locus,
            "is_selected": index == TEST_LOCUS_INDEX,
            "akbari_z": akbari["akbari_z"],
            "akbari_p": akbari["akbari_p"],
            "akbari_converged": akbari["akbari_converged"],
            "akbari_score": abs(akbari["akbari_z"]) if np.isfinite(akbari["akbari_z"]) else np.nan,
            "traj_all_z": float(row.z_score),
            "traj_all_p": float(row.p_value),
            "traj_all_score": traj_score,
            "traj_max_subset_score": best,
            "traj_max_subset_size": best_size,
        })
    frame = pd.DataFrame.from_records(records)
    frame["replicate"] = int(seed)
    return frame


# ------------------------------------------------------------------------- metrics
def _threshold(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    return float(np.quantile(values, 1.0 - ALPHA))


def _score_column(method: str) -> str:
    return {"akbari": "akbari_score", "traj_all": "traj_all_score",
            "traj_max_subset": "traj_max_subset_score"}[method]


def _p_column(method: str) -> str | None:
    return {"akbari": "akbari_p", "traj_all": "traj_all_p",
            "traj_max_subset": None}[method]


def calibration_thresholds(neutral_frame: pd.DataFrame) -> dict[str, float]:
    """95th-percentile neutral-locus scores for each statistic (matched null)."""
    neutral = neutral_frame[~neutral_frame.is_selected]
    thresholds = {}
    for method in METHODS:
        column = _score_column(method)
        rows = neutral
        if method == "akbari":
            rows = neutral[neutral.akbari_converged == 1.0]
        thresholds[method] = _threshold(rows[column].to_numpy(dtype=float))
    return thresholds


def _cell_metric_rows(cell: pd.DataFrame, scenario: str, f: float, strength: float,
                      n_selected: int, thresholds: dict[str, float]) -> list[dict[str, Any]]:
    neutral = cell[~cell.is_selected]
    selected = cell[cell.is_selected]
    naive_device_threshold = thresholds.get("traj_all", float("nan"))
    rows = []
    for method in METHODS:
        score_col = _score_column(method)
        p_col = _p_column(method)
        score_threshold = thresholds.get(method, float("nan"))

        def _subset(frame: pd.DataFrame) -> pd.DataFrame:
            if method == "akbari":
                return frame[frame.akbari_converged == 1.0]
            return frame

        base = _subset(neutral)
        fpr_nominal = fpr_calibrated = np.nan
        if p_col is not None:
            p = base[p_col].to_numpy(dtype=float)
            p = p[np.isfinite(p)]
            if p.size:
                fpr_nominal = float((p < ALPHA).mean())
        scores = base[score_col].to_numpy(dtype=float)
        scores = scores[np.isfinite(scores)]
        n_tests = int(scores.size)
        if scores.size and np.isfinite(score_threshold):
            fpr_calibrated = float((scores > score_threshold).mean())
        fpr_naive_device = np.nan
        if method == "traj_max_subset" and scores.size and np.isfinite(naive_device_threshold):
            fpr_naive_device = float((scores > naive_device_threshold).mean())

        power_nominal = power_calibrated = power_naive_device = np.nan
        if n_selected > 0 and len(selected):
            sel = _subset(selected)
            if p_col is not None:
                p = sel[p_col].to_numpy(dtype=float)
                p = p[np.isfinite(p)]
                if p.size:
                    power_nominal = float((p < ALPHA).mean())
            sel_scores = sel[score_col].to_numpy(dtype=float)
            sel_scores = sel_scores[np.isfinite(sel_scores)]
            if sel_scores.size and np.isfinite(score_threshold):
                power_calibrated = float((sel_scores > score_threshold).mean())
            if (method == "traj_max_subset" and sel_scores.size
                    and np.isfinite(naive_device_threshold)):
                power_naive_device = float((sel_scores > naive_device_threshold).mean())

        if method == "akbari":
            finite = cell.akbari_converged.to_numpy(dtype=float)
            convergence = float(np.nanmean(finite)) if np.isfinite(finite).any() else np.nan
        else:
            convergence = float(np.isfinite(cell[score_col].to_numpy(dtype=float)).mean())
        subset_sizes = cell["traj_max_subset_size"].to_numpy(dtype=float)
        mean_best_size = float(np.mean(subset_sizes)) if method == "traj_max_subset" else np.nan
        rows.append({
            "scenario": scenario, "f": f, "strength": strength,
            "n_trajectories": (int(cell["n_trajectories"].iloc[0])
                               if "n_trajectories" in cell.columns else np.nan),
            "n_selected": n_selected, "method": method, "n_tests": n_tests,
            "fpr_nominal": fpr_nominal, "fpr_calibrated": fpr_calibrated,
            "fpr_naive_device": fpr_naive_device,
            "power_nominal": power_nominal, "power_calibrated": power_calibrated,
            "power_naive_device": power_naive_device,
            "mean_best_subset_size": mean_best_size, "convergence_rate": convergence,
        })
    return rows


# ---------------------------------------------------------------------------- driver
def run_sweep(*, runs: int, null_runs: int, n_loci: int, samples: int,
              akbari_loci: int, f_levels: Sequence[float], strength_levels: Sequence[float],
              time_budget: float | None = None) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run the (f, strength) grid plus the two FPR regimes; write CSV incrementally."""
    base = load_config(SCENARIO_PATH)
    k = len(base["populations"])
    scenario_seed = int(base["seed"])
    started = time.time()
    metric_rows: list[dict[str, Any]] = []
    design: dict[str, Any] = {
        "scenario": str(SCENARIO_PATH.name),
        "n_trajectories": k,
        "runs": int(runs),
        "null_runs": int(null_runs),
        "n_loci": int(n_loci),
        "samples_per_population": int(samples),
        "akbari_loci": int(min(akbari_loci, n_loci)),
        "alpha": ALPHA,
        "baseline_effect": BASELINE_EFFECT,
        "f_levels": [float(f) for f in f_levels],
        "strength_levels": [float(s) for s in strength_levels],
        "min_subset_size": MIN_SUBSET_SIZE,
        "subset_device_definition": "max over subsets |S|>=2 of |meta z| * sign concordance",
        "replicate_seed_formula": "scenario_seed * 10000 + replicate",
        "calibration": "95th percentile neutral-locus score, K=6 no-selection regime",
        "null_regimes": {
            "neutral_drift": "K=6, no selection, no replacement",
            "replacement_only": "K=6, no selection, replacement pop1<-pop2 g12 f=0.85",
        },
        "partial": False,
    }

    def _collect(config: dict[str, Any], scenario: str, f: float, strength: float,
                 n_selected: int, seeds: Iterable[int]) -> pd.DataFrame:
        frames = []
        for replicate, seed in enumerate(seeds):
            frame = run_replicate(config, seed, akbari_loci)
            frame["scenario"] = scenario
            frame["f"] = f
            frame["strength"] = strength
            frame["n_trajectories"] = k
            frames.append(frame)
        return pd.concat(frames, ignore_index=True)

    # --- matched neutral reference (calibration) --------------------------------
    neutral_config = subset_config(base, n_loci, samples, 0, 0.0)
    neutral_seeds = [scenario_seed * 10_000 + r for r in range(null_runs)]
    neutral_frame = _collect(neutral_config, "neutral_drift", np.nan, np.nan, 0, neutral_seeds)
    thresholds = calibration_thresholds(neutral_frame)
    design["calibration_thresholds"] = thresholds
    metric_rows.extend(_cell_metric_rows(neutral_frame, "neutral_drift", np.nan, np.nan,
                                         0, thresholds))

    # --- replacement_only regime (FPR cost) -------------------------------------
    if not (time_budget is not None and (time.time() - started) > time_budget):
        replacement_config = subset_config(base, n_loci, samples, 0, 0.0, replace=True)
        replacement_frame = _collect(replacement_config, "replacement_only", np.nan, np.nan,
                                     0, neutral_seeds)
        metric_rows.extend(_cell_metric_rows(replacement_frame, "replacement_only", np.nan,
                                             np.nan, 0, thresholds))

    # --- (f, strength) selection grid -------------------------------------------
    for f in f_levels:
        n_selected = max(1, _fraction_to_count(f, k))
        for strength in strength_levels:
            if time_budget is not None and (time.time() - started) > time_budget:
                design["partial"] = True
                design["stop_reason"] = "time budget reached"
                break
            config = subset_config(base, n_loci, samples, n_selected, strength)
            seeds = [scenario_seed * 10_000 + r for r in range(runs)]
            cell = _collect(config, "subset_selection", float(f), float(strength),
                            n_selected, seeds)
            metric_rows.extend(_cell_metric_rows(cell, "subset_selection", float(f),
                                                 float(strength), n_selected, thresholds))
            design["wall_seconds"] = round(time.time() - started, 2)
            _atomic_write_csv(pd.DataFrame(metric_rows, columns=list(CSV_COLUMNS)), CSV_PATH)
            _atomic_write_text(SUMMARY_PATH,
                               json.dumps(dict(design, n_rows=len(metric_rows)),
                                          indent=2, sort_keys=True, default=str) + "\n")
        if design.get("partial"):
            break

    metrics = pd.DataFrame(metric_rows, columns=list(CSV_COLUMNS))
    design["wall_seconds"] = round(time.time() - started, 2)
    _atomic_write_csv(metrics, CSV_PATH)
    design["claims"] = _claims(metrics)
    _atomic_write_text(SUMMARY_PATH,
                       json.dumps(dict(design, n_rows=len(metrics)),
                                  indent=2, sort_keys=True, default=str) + "\n")
    return metrics, design


def _claims(metrics: pd.DataFrame) -> list[str]:
    """Deterministic claims read straight off the aggregated CSV."""
    claims: list[str] = []
    grid = metrics[metrics.scenario == "subset_selection"]
    if grid.empty:
        return ["observation: no subset_selection cells were scored."]
    winners = []
    for (f, strength), cell in grid.groupby(["f", "strength"], sort=True):
        by_method = cell.set_index("method")
        device = by_method.loc["traj_max_subset"]
        akbari = by_method.loc["akbari"]
        if (device.power_calibrated > 0.5) and (akbari.power_calibrated <= 0.5):
            winners.append((f, strength, float(device.power_calibrated),
                            float(akbari.power_calibrated)))
    if winners:
        text = "; ".join(
            f"f={f:g}, strength={s:g} (device {d:.2f} vs akbari {a:.2f})"
            for f, s, d, a in winners)
        claims.append("observation: max-over-subsets beats the GLMM (power>0.5 vs <=0.5) at " + text)
    else:
        claims.append("observation: no (f, strength) cell has device power>0.5 while akbari<=0.5.")
    for regime in ("neutral_drift", "replacement_only"):
        rows = metrics[(metrics.scenario == regime) & (metrics.method == "traj_max_subset")]
        if rows.empty:
            continue
        row = rows.iloc[0]
        claims.append(
            f"observation: {regime} max-over-subsets FPR is {row.fpr_calibrated:.3f} "
            f"calibrated vs {row.fpr_naive_device:.3f} with the single-threshold "
            f"(no multiplicity correction) device.")
    return claims


def make_figure(metrics: pd.DataFrame, figure_path: Path) -> None:
    """Power-vs-f lines per strength and FPR bars for the null regimes."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    grid = metrics[metrics.scenario == "subset_selection"]
    strengths = sorted(grid.strength.dropna().unique())
    method_style = {"akbari": ("o", "akbari GLMM"), "traj_all": ("s", "traj all"),
                    "traj_max_subset": ("^", "traj max-subset")}
    fig, axes = plt.subplots(1, max(1, len(strengths)) + 1, figsize=(4 * (len(strengths) + 1), 4.2))
    if len(strengths) == 0:
        strengths = [1.0]
    for ax, strength in zip(axes, strengths):
        cell = grid[grid.strength == strength]
        for method, (marker, label) in method_style.items():
            rows = cell[cell.method == method].sort_values("f")
            ax.plot(rows.f, rows.power_calibrated, marker=marker, label=label)
        ax.axhline(0.5, color="grey", lw=0.8, ls="--")
        ax.set_xlabel("fraction of trajectories selected")
        ax.set_ylabel("calibrated power")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title(f"strength = {strength:g}")
        ax.legend(fontsize=7)
    ax = axes[-1]
    null_rows = metrics[metrics.scenario.isin(("neutral_drift", "replacement_only"))]
    labels, values, colors = [], [], []
    for regime in ("neutral_drift", "replacement_only"):
        for method in ("akbari", "traj_all", "traj_max_subset"):
            row = null_rows[(null_rows.scenario == regime) & (null_rows.method == method)]
            if row.empty:
                continue
            value = row.iloc[0].fpr_calibrated
            if np.isfinite(value):
                labels.append(f"{regime}\n{method}")
                values.append(value)
                colors.append("#333333" if method != "traj_max_subset" else "#b30000")
    ax.bar(range(len(values)), values, color=colors)
    ax.axhline(ALPHA, color="grey", lw=0.8, ls="--")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=6)
    ax.set_ylabel("calibrated FPR")
    ax.set_title("Null-regime FPR")
    fig.tight_layout()
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_path, dpi=150)
    plt.close(fig)


QUICK_DEFAULTS = {"runs": 6, "null_runs": 8, "n_loci": 10, "samples": 6,
                  "akbari_loci": 3, "f_levels": (0.2, 0.4, 0.6, 1.0), "strength_levels": (1.0,)}
FULL_DEFAULTS = {"runs": 20, "null_runs": 40, "n_loci": 20, "samples": 8,
                 "akbari_loci": 6, "f_levels": F_LEVELS, "strength_levels": STRENGTH_LEVELS}


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--quick", action="store_true", help="bounded validation run (<60s)")
    mode.add_argument("--full", action="store_true", help="report run")
    parser.add_argument("--runs", type=int, default=None)
    parser.add_argument("--null-runs", type=int, default=None)
    parser.add_argument("--n-loci", type=int, default=None)
    parser.add_argument("--samples-per-population", type=int, default=None)
    parser.add_argument("--akbari-loci", type=int, default=None)
    parser.add_argument("--time-budget", type=float, default=None)
    parser.add_argument("--output", type=Path, default=CSV_PATH)
    parser.add_argument("--figure", type=Path, default=FIGURE_PATH)
    parser.add_argument("--no-figure", action="store_true")
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    quick = bool(args.quick)
    defaults = QUICK_DEFAULTS if quick else FULL_DEFAULTS
    global CSV_PATH, SUMMARY_PATH
    if args.output != CSV_PATH:
        CSV_PATH = args.output
        SUMMARY_PATH = args.output.with_name("summary.json")
    runs = args.runs if args.runs is not None else defaults["runs"]
    null_runs = args.null_runs if args.null_runs is not None else defaults["null_runs"]
    n_loci = args.n_loci if args.n_loci is not None else defaults["n_loci"]
    samples = (args.samples_per_population if args.samples_per_population is not None
               else defaults["samples"])
    akbari_loci = args.akbari_loci if args.akbari_loci is not None else defaults["akbari_loci"]
    budget = args.time_budget if args.time_budget is not None else (55.0 if quick else 660.0)
    metrics, design = run_sweep(
        runs=runs, null_runs=null_runs, n_loci=n_loci, samples=samples,
        akbari_loci=akbari_loci, f_levels=defaults["f_levels"],
        strength_levels=defaults["strength_levels"], time_budget=budget)
    label = "quick" if quick else "full"
    print(f"[{label}] {len(metrics)} metric rows in {design['wall_seconds']}s; "
          f"partial={design['partial']}")
    print(f"thresholds: {json.dumps(design.get('calibration_thresholds', {}))}")
    for claim in design.get("claims", []):
        print(claim)
    if not args.no_figure and not metrics.empty:
        make_figure(metrics, args.figure)
        print(f"figure: {args.figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
