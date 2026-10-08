#!/usr/bin/env python3
"""Seeded multi-stratum neutral holdout check of empirical matched-control ranks.

This is a binomial/shared-drift null stress test, not an ancient-DNA demographic
validation or proof that QC metadata are correctly estimated in real cohorts.
"""

import argparse
import json

import numpy as np
import pandas as pd

from scripts.provenance import RunRecorder
from sim.empirical_null import empirical_null


SEED_HELP = "explicit PCG64 seed; no global or unseeded draws"


def simulate_null(seed: int) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Independent neutral loci across four differentiated genomic strata."""
    rng = np.random.default_rng(seed)
    counts, metadata, candidates = [], [], []
    for stratum in range(4):
        chrom = str(1 + stratum // 2)
        high = stratum % 2
        base = 0.24 if high == 0 else 0.68
        times = (0, 1, 2, 3) if high == 0 else (1, 2, 3, 4)
        for index in range(80):  # 48 controls and 32 prespecified held-out neutral loci
            locus = f"s{stratum}_{index:02d}"
            if index >= 48:
                candidates.append(locus)
            metadata.append(dict(locus=locus, chromosome=chrom,
                                 position=(80 * stratum + index + 1) * 2_000_000,
                                 allele_frequency=base + rng.uniform(-.004, .004),
                                 recombination_rate=(.8 if high == 0 else 2.2) + rng.uniform(-.02, .02),
                                 local_ld=(.18 if high == 0 else .65) + rng.uniform(-.005, .005),
                                 imputation_quality=.96 + rng.uniform(-.004, .004),
                                 missingness=.025 + rng.uniform(-.004, .004),
                                 temporal_coverage=(.85 if high == 0 else .95)
                                 + rng.uniform(-.004, .004)))
            # Shared neutral time shocks induce correlated trajectories within
            # a locus; they are generated independently across loci.
            shared = rng.normal(0, .018, size=4).cumsum()
            origin = base + rng.normal(0, .015)
            for trajectory, offset in (("a", -.05), ("b", 0), ("c", .05)):
                local = rng.normal(0, .016, size=4).cumsum()
                for j, time in enumerate(times):
                    depth = 70 if high == 0 else 130
                    probability = float(np.clip(origin + offset + shared[j] + local[j], .02, .98))
                    counts.append((locus, trajectory, time,
                                   int(rng.binomial(depth, probability)), depth))
    return (pd.DataFrame(counts, columns=["locus", "trajectory", "time", "alt_count", "total_count"]),
            pd.DataFrame(metadata), candidates)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True, help=SEED_HELP)
    args = parser.parse_args()
    with RunRecorder("sim-ctrl-null-validation", command=[
        "python3", "-m", "sim.validate_empirical_null", "--seed", str(args.seed)
    ]) as run:
        run.set_seed(args.seed)
        for source in ("sim/empirical_null.py", "sim/trajectory_stat.py",
                       "sim/validate_empirical_null.py"):
            run.record_input(source)
        run.record_parameter("design", "4 strata x (48 neutral controls + 32 neutral candidates); "
                             "3 trajectories, 4 times; independent loci, correlated neutral drift within locus")
        run.record_parameter("matcher", "default calipers, same chromosome and exact time support; "
                             "LD exclusion 1Mb; min_controls=20, max_controls=40")
        counts, metadata, candidates = simulate_null(args.seed)
        for name, data in (("counts", counts), ("metadata", metadata)):
            path = run.run_dir / f"{name}.csv"
            data.to_csv(path, index=False)
            run.record_output(path)
        null = empirical_null(counts, metadata, candidates, min_controls=20, max_controls=40)
        for name, data in (("ranks", null.results), ("matches", null.matches)):
            path = run.run_dir / f"{name}.csv"
            data.to_csv(path, index=False)
            run.record_output(path)
        p = np.sort(null.results.p_empirical.to_numpy())
        n = len(p)
        ks = float(max(np.max(np.arange(1, n + 1) / n - p),
                       np.max(p - np.arange(n) / n)))
        metrics = {"seed": args.seed, "n_candidates": n, "controls_per_candidate": 40,
                   "mean_p": float(p.mean()), "ks_uniform_d": ks,
                   "fraction_p_le_0_05": float((p <= .05).mean()),
                   "fraction_p_le_0_10": float((p <= .10).mean())}
        path = run.run_dir / "calibration.json"
        path.write_text(json.dumps(metrics, indent=2) + "\n")
        run.record_output(path)
        print(json.dumps(metrics, indent=2))
        assert n == 128 and (null.results.n_controls == 40).all()
        assert .43 <= metrics["mean_p"] <= .58
        assert metrics["ks_uniform_d"] < .15
        assert metrics["fraction_p_le_0_05"] <= .10
        assert metrics["fraction_p_le_0_10"] <= .16


if __name__ == "__main__":
    main()
