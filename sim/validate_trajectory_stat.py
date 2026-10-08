#!/usr/bin/env python3
"""Run the fixed-design sign recovery check and record reproducible outputs."""

import argparse

import numpy as np
import pandas as pd

from scripts.provenance import RunRecorder
from sim.trajectory_stat import trajectory_statistic


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args()
    with RunRecorder("sim-stat-validation", command=[
        "python3", "-m", "sim.validate_trajectory_stat", "--seed", str(args.seed)
    ]) as run:
        run.set_seed(args.seed)
        for source in ("sim/trajectory_stat.py", "sim/validate_trajectory_stat.py"):
            run.record_input(source)
        run.record_parameter("design", "three trajectories, nine times, depth 2000; slopes +/-0.45")
        rng = np.random.default_rng(args.seed)
        rows = []
        for locus, trend in (("up", 0.45), ("down", -0.45)):
            for trajectory, intercept in (("a", 0.3), ("b", 0.45), ("c", 0.6)):
                for time in np.linspace(0, 1, 9):
                    probability = intercept + trend * (time - 0.5)
                    rows.append((locus, trajectory, time,
                                 rng.binomial(2000, probability), 2000))
        data = pd.DataFrame(rows, columns=("locus", "trajectory", "time",
                                           "alt_count", "total_count"))
        result = trajectory_statistic(data)
        assert (result.slopes.query("locus == 'up'").slope > 0.3).all()
        assert (result.slopes.query("locus == 'down'").slope < -0.3).all()
        assert set(("locus", "trajectory", "slope", "se", "direction", "n", "span")) <= set(result.slopes)
        for name, frame in (("counts", data), ("slopes", result.slopes),
                            ("summary", result.summary)):
            path = run.run_dir / f"{name}.csv"
            frame.to_csv(path, index=False)
            run.record_output(path)
        print(result.summary[["locus", "concordance", "meta_slope", "meta_se"]].to_string(index=False))


if __name__ == "__main__":
    main()
