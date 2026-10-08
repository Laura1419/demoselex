#!/usr/bin/env python3
"""
Run all simulation scenarios and record provenance
"""

import os
import sys
import json
import yaml
from pathlib import Path
from simulator import DemographicSimulator, load_config, save_results
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../..'))
from scripts.provenance import Provenance

def run_all_scenarios(config_file: str, output_base_dir: str):
    """Run all scenarios defined in the configuration"""

    # Load configuration
    config = load_config(config_file)
    base_seed = config.get('seed', 42)

    # Create output directory
    output_base_dir = Path(output_base_dir)
    output_base_dir.mkdir(parents=True, exist_ok=True)

    # Run each scenario
    scenarios = config.get('scenarios', {})
    all_results = {}

    for scenario_name, scenario_config in scenarios.items():
        print(f"\n=== Running scenario: {scenario_name} ===")

        # Set up scenario-specific output directory
        scenario_output_dir = output_base_dir / scenario_name

        # Set seed for reproducibility
        scenario_config['seed'] = base_seed + hash(scenario_name) % 1000

        # Record provenance for this run
        prov = Provenance(
            seed=scenario_config['seed'],
            config_file=config_file,
            output_dir=str(scenario_output_dir)
        )

        # Run simulation
        simulator = DemographicSimulator(config)
        simulator.setup_populations()
        results = simulator.run_scenario(scenario_config)

        # Save results
        save_results(results, scenario_output_dir)

        # Record provenance
        output_files = {
            'config': config_file,
            'genotype_matrix_pop1': str(scenario_output_dir / 'pop1_genotypes.csv'),
            'genotype_matrix_pop2': str(scenario_output_dir / 'pop2_genotypes.csv'),
            'truth_table': str(scenario_output_dir / 'truth_table.json'),
            'summary': str(scenario_output_dir / 'summary.json')
        }

        provenance_file = prov.record_run(output_files, scenario_config)
        print(f"Results saved to: {scenario_output_dir}")
        print(f"Provenance recorded to: {provenance_file}")

        all_results[scenario_name] = {
            'output_dir': str(scenario_output_dir),
            'provenance_file': provenance_file
        }

    return all_results

def create_scenario_summary(results: dict, output_dir: str):
    """Create a summary of all scenarios run"""
    summary = {
        'scenarios_run': list(results.keys()),
        'total_scenarios': len(results),
        'scenario_details': {}
    }

    for scenario_name, scenario_info in results.items():
        # Load summary for each scenario
        summary_file = Path(scenario_info['output_dir']) / 'summary.json'
        if summary_file.exists():
            with open(summary_file, 'r') as f:
                scenario_summary = json.load(f)
                summary['scenario_details'][scenario_name] = scenario_summary

    # Save overall summary
    summary_file = Path(output_dir) / 'simulation_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\n=== Overall Summary ===")
    print(f"Total scenarios run: {summary['total_scenarios']}")
    print(f"Scenarios: {', '.join(summary['scenarios_run'])}")
    print(f"Summary saved to: {summary_file}")

def main():
    """Main function to run simulations"""
    if len(sys.argv) != 3:
        print("Usage: python3 run_simulations.py <config_file> <output_dir>")
        sys.exit(1)

    config_file = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(config_file):
        print(f"Config file not found: {config_file}")
        sys.exit(1)

    print("Starting demographic simulations...")
    print(f"Config: {config_file}")
    print(f"Output: {output_dir}")

    # Run all scenarios
    results = run_all_scenarios(config_file, output_dir)

    # Create summary
    create_scenario_summary(results, output_dir)

    print("\n=== Simulation Complete ===")
    print(f"All scenarios completed successfully!")

if __name__ == "__main__":
    main()