#!/usr/bin/env python3
"""Run one or every fixed-seed scenario and record the repository provenance manifest.

Usage: python3 sim/scripts/run_simulations.py sim/scenarios sim/output/demo
       python3 sim/scripts/run_simulations.py sim/scenarios/neutral_drift.yaml sim/output/demo
"""

from __future__ import annotations

import argparse
from pathlib import Path
import shutil
import sys

# Support invocation from the repository root without a package install.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from scripts.provenance import RunRecorder  # noqa: E402
from simulator import DemographicSimulator, load_config, save_results  # noqa: E402


def run(config_path: Path, output_base: Path) -> Path:
    config = load_config(config_path)
    if config_path.stem != config.get("type"):
        raise ValueError(f"Config name/type mismatch: {config_path}")
    simulator = DemographicSimulator(config)
    destination = output_base / config["type"]
    with RunRecorder(f"sim-{config['type']}-{config['seed']}", command=sys.argv) as recorder:
        recorder.set_seed(config["seed"])
        recorder.record_input(config_path)
        for key in ("type", "generations", "n_loci", "populations", "migration",
                    "replacement", "selection", "recombination_rate", "sample_generations"):
            recorder.record_parameter(key, config.get(key))
        paths = save_results(simulator.run(), destination)
        for path in paths:
            recorder.record_output(path)
    # Copy the finalized manifest alongside the matrices for downstream consumers.
    manifest = destination / "provenance.json"
    shutil.copyfile(recorder.run_dir / "manifest.json", manifest)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path, help="A scenario .yaml or directory of scenarios")
    parser.add_argument("output", type=Path, help="Output parent directory")
    args = parser.parse_args()
    paths = sorted(args.config.glob("*.yaml")) if args.config.is_dir() else [args.config]
    if not paths:
        parser.error("No scenario YAML files found")
    for path in paths:
        print(f"{path.stem}: {run(path, args.output)}")


if __name__ == "__main__":
    main()
