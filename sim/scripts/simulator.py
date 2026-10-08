"""Seeded, forward-time diploid mosaic/Wright–Fisher demographic simulator.

This is a small validation simulator, not an inference model or SLiM substitute:
parents are drawn with additive viability weights, independent gametes are
recombined along an ordered panel, and whole diploid individuals migrate.
All populations advance synchronously; samples are snapshots, not ancestors.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml


SCENARIOS = {
    "neutral_drift", "replacement_only", "one_population_selection",
    "shared_selection", "opposing_selection", "selection_aligned_turnover",
}


def load_config(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, dict):
        raise ValueError("Scenario config must be a mapping")
    return config


class DemographicSimulator:
    """Simulate one scenario; constructing/running a new instance resets its RNG."""

    def __init__(self, config: dict[str, Any]):
        self.config = config
        self._validate()
        self.seed = config["seed"]
        self.rng = np.random.default_rng(self.seed)
        self.n_loci = config["n_loci"]
        self.generations = config["generations"]
        self.names = list(config["populations"])
        self.populations: dict[str, np.ndarray] = {}

    def _validate(self) -> None:
        c = self.config
        if c.get("type") not in SCENARIOS:
            raise ValueError(f"Unknown scenario type: {c.get('type')}")
        if type(c.get("seed")) is not int or c["seed"] < 0:
            raise ValueError("An explicit nonnegative integer seed is required")
        if not isinstance(c.get("n_loci"), int) or c["n_loci"] < 2:
            raise ValueError("n_loci must be at least 2")
        if not isinstance(c.get("generations"), int) or c["generations"] < 1:
            raise ValueError("generations must be positive")
        if not isinstance(c.get("populations"), dict) or len(c["populations"]) < 2:
            raise ValueError("At least two named populations are required")
        for name, pop in c["populations"].items():
            if not isinstance(pop.get("size"), int) or pop["size"] < 2:
                raise ValueError(f"{name}: size must be at least 2")
        for key in ("recombination_rate",):
            if not 0 <= c.get(key, 0.02) <= 1:
                raise ValueError(f"{key} must be between 0 and 1")
        if c.get("generation_years", 25) <= 0:
            raise ValueError("generation_years must be positive")
        if not isinstance(c.get("sample_generations"), list) or not c["sample_generations"]:
            raise ValueError("sample_generations must be a nonempty list")
        if any(type(g) is not int or g < 0 or g > c["generations"] for g in c["sample_generations"]):
            raise ValueError("sample_generations must fall within 0..generations")
        if len(set(c["sample_generations"])) != len(c["sample_generations"]):
            raise ValueError("sample_generations must be distinct")
        if not isinstance(c.get("samples_per_population"), int) or c["samples_per_population"] < 1:
            raise ValueError("samples_per_population must be positive")
        names = set(c["populations"])
        for pop, loci in c.get("selection", {}).items():
            if pop not in names:
                raise ValueError(f"Unknown selected population: {pop}")
            for locus, effect in loci.items():
                if not 0 <= int(locus) < c["n_loci"] or effect <= -0.5:
                    raise ValueError("Invalid locus or additive selection effect")
        for recipient, donors in c.get("migration", {}).items():
            if recipient not in names or any(d not in names or d == recipient for d in donors):
                raise ValueError("Migration donor/recipient must be distinct named populations")
            if any(rate < 0 or rate > 1 for rate in donors.values()) or sum(donors.values()) > 1:
                raise ValueError("Migration rates must be nonnegative and sum to <= 1")
        for event in c.get("replacement", []):
            if (event["recipient"] not in names or event["source"] not in names
                    or event["recipient"] == event["source"]
                    or not 1 <= event["generation"] <= c["generations"]
                    or not 0 <= event["fraction"] <= 1):
                raise ValueError("Invalid replacement pulse")
        for pop, loci in c.get("founder_frequencies", {}).items():
            if pop not in names or any(not 0 <= int(k) < c["n_loci"] or not 0 <= v <= 1
                                       for k, v in loci.items()):
                raise ValueError("Invalid founder frequency")

    def _initialize(self) -> None:
        # Shared ancestral frequencies plus population-specific divergence.
        ancestral = self.rng.beta(2, 2, self.n_loci)
        for name, pop in self.config["populations"].items():
            frequencies = self.rng.beta(8 * ancestral, 8 * (1 - ancestral))
            for locus, value in self.config.get("founder_frequencies", {}).get(name, {}).items():
                frequencies[int(locus)] = value
            self.populations[name] = (self.rng.random((pop["size"], 2, self.n_loci))
                                      < frequencies).astype(np.uint8)

    def _gametes(self, individuals: np.ndarray, weights: np.ndarray, count: int) -> np.ndarray:
        """Sample parents with replacement and recombine their two homologues."""
        parent = self.rng.choice(len(individuals), size=count, p=weights)
        start = self.rng.integers(0, 2, size=(count, 1), dtype=np.int8)
        switches = self.rng.random((count, self.n_loci - 1)) < self.config.get("recombination_rate", 0.02)
        strands = np.concatenate((start, (start + np.cumsum(switches, axis=1)) % 2), axis=1)
        return individuals[parent[:, None], strands, np.arange(self.n_loci)[None, :]]

    def _reproduce(self, name: str, individuals: np.ndarray) -> np.ndarray:
        effects = self.config.get("selection", {}).get(name, {})
        fitness = np.ones(len(individuals), dtype=float)
        for locus, effect in effects.items():
            fitness += effect * individuals[:, :, int(locus)].sum(axis=1)
        # effects > -0.5 per allele, hence fitness stays strictly positive for one locus;
        # for multiple negative loci a floor prevents invalid probabilities.
        fitness = np.maximum(fitness, 1e-9)
        fitness /= fitness.sum()
        gametes = self._gametes(individuals, fitness, 2 * len(individuals))
        return gametes.reshape(len(individuals), 2, self.n_loci)

    def _sample(self, generation: int, rows: list[dict[str, Any]]) -> None:
        for name in self.names:
            pool = self.populations[name]
            n = min(len(pool), self.config["samples_per_population"])
            indices = self.rng.choice(len(pool), size=n, replace=False)
            genotypes = pool[indices].sum(axis=1)
            for k, genotype in enumerate(genotypes):
                rows.append({
                    "sample_id": f"{name}_g{generation:04d}_{k:04d}",
                    "trajectory": name,
                    "generation": generation,
                    "years_bp": (self.generations - generation) * self.config.get("generation_years", 25),
                    "genotype": genotype.tolist(),
                })

    def run(self) -> dict[str, Any]:
        self._initialize()
        rows: list[dict[str, Any]] = []
        sample_times = set(self.config["sample_generations"])
        if 0 in sample_times:
            self._sample(0, rows)
        for generation in range(1, self.generations + 1):
            old = self.populations
            next_generation = {name: self._reproduce(name, old[name]) for name in self.names}
            # Migration and pulse replacement draw whole diploid migrants from the
            # previous generation's source, avoiding order dependence among populations.
            for recipient, donors in self.config.get("migration", {}).items():
                pool = next_generation[recipient]
                chosen = self.rng.random(len(pool))
                lower = 0.0
                for source, rate in donors.items():
                    mask = (chosen >= lower) & (chosen < lower + rate)
                    pool[mask] = old[source][self.rng.integers(len(old[source]), size=int(mask.sum()))]
                    lower += rate
            for event in self.config.get("replacement", []):
                if event["generation"] != generation:
                    continue
                pool = next_generation[event["recipient"]]
                mask = self.rng.random(len(pool)) < event["fraction"]
                pool[mask] = old[event["source"]][self.rng.integers(len(old[event["source"]]), size=int(mask.sum()))]
            self.populations = next_generation
            if generation in sample_times:
                self._sample(generation, rows)
        truth = [
            {"locus": locus, "trajectory": name,
             "selected": int(str(locus) in self.config.get("selection", {}).get(name, {})),
             "effect_size": self.config.get("selection", {}).get(name, {}).get(str(locus), 0.0)}
            for name in self.names for locus in range(self.n_loci)
        ]
        return {"samples": rows, "truth": truth, "type": self.config["type"], "seed": self.seed}


def save_results(results: dict[str, Any], output_dir: str | Path) -> list[Path]:
    """Write an aligned individual-by-locus matrix, metadata, and per-population truth."""
    dest = Path(output_dir)
    dest.mkdir(parents=True, exist_ok=True)
    rows = results["samples"]
    if not rows:
        raise ValueError("No individuals sampled")
    loci = [f"locus_{l}" for l in range(len(rows[0]["genotype"]))]
    matrix, metadata, truth, summary = (dest / name for name in
                                        ("genotypes.csv", "samples.csv", "truth_table.csv", "summary.json"))
    with matrix.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["sample_id", *loci])
        writer.writerows([row["sample_id"], *row["genotype"]] for row in rows)
    with metadata.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["sample_id", "trajectory", "generation", "years_bp"])
        writer.writeheader()
        writer.writerows({key: row[key] for key in writer.fieldnames} for row in rows)
    with truth.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=["locus", "trajectory", "selected", "effect_size"])
        writer.writeheader()
        writer.writerows(results["truth"])
    summary.write_text(json.dumps({"scenario_type": results["type"], "seed": results["seed"],
                                   "n_samples": len(rows), "n_loci": len(loci)}, indent=2) + "\n",
                       encoding="utf-8")
    return [matrix, metadata, truth, summary]
