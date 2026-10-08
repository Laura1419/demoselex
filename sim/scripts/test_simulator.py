"""End-to-end regression of the six declared fixed-seed simulation scenarios."""

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCENARIOS = ROOT / "sim" / "scenarios"
CLI = ROOT / "sim" / "scripts" / "run_simulations.py"
EXPECTED = {
    "neutral_drift": {},
    "replacement_only": {},
    "one_population_selection": {"pop1": 0.25},
    "shared_selection": {"pop1": 0.25, "pop2": 0.25},
    "opposing_selection": {"pop1": 0.25, "pop2": -0.25},
    "selection_aligned_turnover": {"pop2": 0.25},
}


class ScenarioTests(unittest.TestCase):
    def test_all_scenarios_through_cli(self):
        self.assertEqual({p.stem for p in SCENARIOS.glob("*.yaml")}, set(EXPECTED))
        with tempfile.TemporaryDirectory() as temp:
            first, second = Path(temp) / "first", Path(temp) / "second"
            for destination in (first, second):
                subprocess.run([sys.executable, str(CLI), str(SCENARIOS), str(destination)],
                               cwd=ROOT, check=True, capture_output=True, text=True)
            for scenario, effects in EXPECTED.items():
                output = first / scenario
                for name in ("genotypes.csv", "samples.csv", "truth_table.csv", "summary.json"):
                    self.assertEqual((output / name).read_bytes(), (second / scenario / name).read_bytes())
                with (output / "genotypes.csv").open() as f:
                    genotypes = list(csv.DictReader(f))
                with (output / "samples.csv").open() as f:
                    samples = list(csv.DictReader(f))
                with (output / "truth_table.csv").open() as f:
                    truth = list(csv.DictReader(f))
                self.assertEqual(len(genotypes), len(samples))
                self.assertEqual(len(genotypes), 160)
                self.assertEqual([r["sample_id"] for r in genotypes],
                                 [r["sample_id"] for r in samples])
                self.assertEqual(set(r["years_bp"] for r in samples), {"0", "200", "400", "600"})
                self.assertEqual(len(truth), 256)
                self.assertTrue(all(0 <= int(row[f"locus_{i}"]) <= 2
                                    for row in genotypes for i in range(128)))
                selected = {row["trajectory"]: float(row["effect_size"]) for row in truth
                            if row["locus"] == "0" and row["selected"] == "1"}
                self.assertEqual(selected, effects)
                self.assertTrue(all(row["effect_size"] == "0.0" and row["selected"] == "0"
                                    for row in truth if row["locus"] != "0"))
                summary = json.loads((output / "summary.json").read_text())
                self.assertEqual(summary["scenario_type"], scenario)
                manifest = json.loads((output / "provenance.json").read_text())
                self.assertEqual(manifest["seed"], summary["seed"])
                self.assertIn("numpy", manifest["versions"])
                self.assertIn("python", manifest["versions"])
                self.assertIn(str(SCENARIOS / (scenario + ".yaml")), manifest["inputs"])
                self.assertIn(str(output / "genotypes.csv"), manifest["outputs"])
                self.assertIn(str(output / "truth_table.csv"), manifest["outputs"])

    def test_migration_and_replacement_change_lineages(self):
        from simulator import DemographicSimulator, load_config
        config = load_config(SCENARIOS / "replacement_only.yaml")
        result = DemographicSimulator(config).run()
        by_pop_gen = lambda p, g: [x["genotype"][0] for x in result["samples"]
                                  if x["trajectory"] == p and x["generation"] == g]
        frequency = lambda p, g: sum(by_pop_gen(p, g)) / (2 * len(by_pop_gen(p, g)))
        self.assertLess(frequency("pop1", 0), frequency("pop2", 0))
        self.assertGreater(frequency("pop1", 16), frequency("pop1", 0) + 0.3)
        migration = load_config(SCENARIOS / "neutral_drift.yaml")
        migration["migration"] = {"pop1": {"pop2": 1.0}}
        migration["generations"] = 1
        migration["sample_generations"] = [0, 1]
        # Entire recipient generation consists of copies of prior source pool.
        migration["samples_per_population"] = migration["populations"]["pop1"]["size"]
        migrated = DemographicSimulator(migration).run()
        source = {tuple(x["genotype"]) for x in migrated["samples"]
                  if x["trajectory"] == "pop2" and x["generation"] == 0}
        recipient = [tuple(x["genotype"]) for x in migrated["samples"]
                     if x["trajectory"] == "pop1" and x["generation"] == 1]
        self.assertTrue(all(x in source for x in recipient))

    def test_rejects_unseeded_and_invalid_effect(self):
        from simulator import DemographicSimulator, load_config
        config = load_config(SCENARIOS / "neutral_drift.yaml")
        del config["seed"]
        with self.assertRaisesRegex(ValueError, "seed"):
            DemographicSimulator(config)
        config["seed"] = 1
        config["selection"] = {"pop2": {"999": 0.1}}
        with self.assertRaisesRegex(ValueError, "locus"):
            DemographicSimulator(config)


if __name__ == "__main__":
    unittest.main()
