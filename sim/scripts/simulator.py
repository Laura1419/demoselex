#!/usr/bin/env python3
"""
Demographic simulator scaffold
Produces dated individuals with genome-wide genotypes along population trajectories
with configurable migration, replacement, drift, and selection
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
import json
import argparse
from pathlib import Path
import yaml

class Population:
    """Represents a population with individuals, genotypes, and demographic history"""

    def __init__(self, name: str, size: int, n_loci: int, time_origin: int):
        self.name = name
        self.size = size
        self.n_loci = n_loci
        self.time_origin = time_origin

        # Generate initial genotypes (diploid, binary alleles)
        self.individuals = self._initialize_genotypes()
        self.history = []  # Track population history

    def _initialize_genotypes(self) -> np.ndarray:
        """Initialize random genotypes for the population"""
        # Shape: (n_individuals, 2, n_loci) for diploid genotypes
        return np.random.choice([0, 1], size=(self.size, 2, self.n_loci))

    def evolve(self, generations: int, migration_rate: float = 0.0,
               selection_strength: float = 0.0, selected_loci: List[int] = None) -> None:
        """Evolve the population for a number of generations"""
        if selected_loci is None:
            selected_loci = []

        for gen in range(generations):
            # Selection (if specified)
            if selection_strength > 0 and selected_loci:
                self._apply_selection(selection_strength, selected_loci)

            # Drift (reproduction with random mating)
            new_individuals = self._reproduce()

            # Migration (if specified)
            if migration_rate > 0:
                new_individuals = self._apply_migration(new_individuals, migration_rate)

            self.individuals = new_individuals
            self.history.append({
                'generation': gen,
                'size': self.size,
                'average_fitness': self._calculate_average_fitness(selected_loci)
            })

    def _apply_selection(self, strength: float, selected_loci: List[int]) -> None:
        """Apply selection to the population"""
        # Simple additive fitness model
        fitnesses = np.ones(self.size)

        for locus in selected_loci:
            # Count derived alleles (1s) at this locus
            allele_counts = np.sum(self.individuals[:, :, locus], axis=1)
            fitnesses += strength * allele_counts

        # Normalize and use as probabilities for reproduction
        fitnesses = np.maximum(fitnesses, 0)  # Ensure non-negative
        if np.sum(fitnesses) > 0:
            fitnesses /= np.sum(fitnesses)

        # Store for use in reproduction
        self._fitnesses = fitnesses

    def _reproduce(self) -> np.ndarray:
        """Reproduce the population with drift"""
        new_individuals = np.zeros_like(self.individuals)

        # Use fitness-weighted sampling if selection was applied
        if hasattr(self, '_fitnesses') and np.sum(self._fitnesses) > 0:
            parent_indices = np.random.choice(
                self.size, size=self.size, p=self._fitnesses
            )
        else:
            parent_indices = np.random.choice(self.size, size=self.size)

        for i, parent_idx in enumerate(parent_indices):
            # Simple reproduction: choose one chromosome from each parent
            # (for simplicity, using the same parent for both chromosomes)
            new_individuals[i, 0] = self.individuals[parent_idx, np.random.choice(2)]
            new_individuals[i, 1] = self.individuals[parent_idx, np.random.choice(2)]

        return new_individuals

    def _apply_migration(self, individuals: np.ndarray, rate: float) -> np.ndarray:
        """Apply migration by replacing some individuals with migrants"""
        # Simple migration model: replace a fraction of individuals
        n_migrants = int(self.size * rate)
        migrant_indices = np.random.choice(self.size, size=n_migrants, replace=False)

        # Generate random genotypes for migrants
        for idx in migrant_indices:
            individuals[idx] = np.random.choice([0, 1], size=(2, self.n_loci))

        return individuals

    def _calculate_average_fitness(self, selected_loci: List[int]) -> float:
        """Calculate average fitness of the population"""
        if not selected_loci:
            return 1.0

        fitnesses = np.ones(self.size)
        for locus in selected_loci:
            allele_counts = np.sum(self.individuals[:, :, locus], axis=1)
            fitnesses += allele_counts * 0.1  # Fixed selection strength for measurement

        return np.mean(fitnesses)

    def get_genotype_matrix(self) -> np.ndarray:
        """Get genotype matrix (individuals x loci)"""
        # Sum alleles across chromosomes
        return np.sum(self.individuals, axis=1)

    def get_frequencies(self) -> np.ndarray:
        """Get allele frequencies across loci"""
        return np.mean(self.get_genotype_matrix(), axis=0) / 2.0

class DemographicSimulator:
    """Main simulator class that handles multiple populations and scenarios"""

    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.seed = config.get('seed', 42)
        np.random.seed(self.seed)

        self.n_loci = config.get('n_loci', 10000)
        self.populations = {}
        self.trajectories = []

    def setup_populations(self) -> None:
        """Initialize populations based on configuration"""
        pop_configs = self.config.get('populations', {})

        for pop_name, pop_config in pop_configs.items():
            size = pop_config.get('size', 1000)
            time_origin = pop_config.get('time_origin', 0)
            self.populations[pop_name] = Population(pop_name, size, self.n_loci, time_origin)

    def run_scenario(self, scenario_config: Dict[str, Any]) -> Dict[str, Any]:
        """Run a specific scenario and return results"""
        scenario_type = scenario_config['type']

        # Set seed for reproducibility
        np.random.seed(scenario_config.get('seed', self.seed))

        # Reset populations for this scenario
        self.setup_populations()

        if scenario_type == 'neutral_drift':
            return self._run_neutral_drift(scenario_config)
        elif scenario_type == 'replacement_only':
            return self._run_replacement_only(scenario_config)
        elif scenario_type == 'one_population_selection':
            return self._run_one_population_selection(scenario_config)
        elif scenario_type == 'shared_selection':
            return self._run_shared_selection(scenario_config)
        elif scenario_type == 'opposing_selection':
            return self._run_opposing_selection(scenario_config)
        elif scenario_type == 'selection_aligned_turnover':
            return self._run_selection_aligned_turnover(scenario_config)
        else:
            raise ValueError(f"Unknown scenario type: {scenario_type}")

    def _run_neutral_drift(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Neutral drift scenario - no selection, only drift and migration"""
        generations = config.get('generations', 100)
        migration_rate = config.get('migration_rate', 0.01)

        for pop in self.populations.values():
            pop.evolve(generations, migration_rate=migration_rate)

        return self._collect_results(config, selected_loci=[])

    def _run_replacement_only(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Replacement-only scenario - population replacement"""
        # Simple replacement: one population replaces another
        generations = config.get('generations', 100)
        replacement_generation = config.get('replacement_generation', 50)

        for pop_name, pop in self.populations.items():
            if pop_name == 'pop2':
                # This population only exists after replacement
                if generations > replacement_generation:
                    remaining_generations = generations - replacement_generation
                    pop.evolve(remaining_generations)
            else:
                pop.evolve(generations)

        return self._collect_results(config, selected_loci=[])

    def _run_one_population_selection(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Selection in one population only"""
        generations = config.get('generations', 100)
        selection_strength = config.get('selection_strength', 0.1)
        selected_loci = config.get('selected_loci', [0])  # Default: locus 0

        for pop_name, pop in self.populations.items():
            if pop_name == config.get('selected_population', 'pop1'):
                pop.evolve(generations, selection_strength=selection_strength,
                          selected_loci=selected_loci)
            else:
                pop.evolve(generations)

        return self._collect_results(config, selected_loci)

    def _run_shared_selection(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Shared selection across populations"""
        generations = config.get('generations', 100)
        selection_strength = config.get('selection_strength', 0.1)
        selected_loci = config.get('selected_loci', [0])

        for pop in self.populations.values():
            pop.evolve(generations, selection_strength=selection_strength,
                      selected_loci=selected_loci)

        return self._collect_results(config, selected_loci)

    def _run_opposing_selection(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Opposing selection in different populations"""
        generations = config.get('generations', 100)
        selection_strength = config.get('selection_strength', 0.1)
        selected_loci = config.get('selected_loci', [0])

        for pop_name, pop in self.populations.items():
            if pop_name == 'pop1':
                # Positive selection in pop1
                pop.evolve(generations, selection_strength=selection_strength,
                          selected_loci=selected_loci)
            elif pop_name == 'pop2':
                # Negative selection in pop2 (opposing direction)
                pop.evolve(generations, selection_strength=-selection_strength,
                          selected_loci=selected_loci)
            else:
                pop.evolve(generations)

        return self._collect_results(config, selected_loci)

    def _run_selection_aligned_turnover(self, config: Dict[str, Any]) -> Dict[str, Any]:
        """Selection aligned with population turnover"""
        generations = config.get('generations', 100)
        selection_strength = config.get('selection_strength', 0.1)
        selected_loci = config.get('selected_loci', [0])
        replacement_generation = config.get('replacement_generation', 50)

        for pop_name, pop in self.populations.items():
            if pop_name == 'pop2':
                # Population with aligned selection
                if generations > replacement_generation:
                    remaining_generations = generations - replacement_generation
                    pop.evolve(remaining_generations, selection_strength=selection_strength,
                              selected_loci=selected_loci)
            else:
                pop.evolve(generations)

        return self._collect_results(config, selected_loci)

    def _collect_results(self, config: Dict[str, Any], selected_loci: List[int]) -> Dict[str, Any]:
        """Collect results from all populations"""
        results = {
            'scenario_type': config['type'],
            'seed': config.get('seed', self.seed),
            'selected_loci': selected_loci,
            'populations': {},
            'truth_table': {}
        }

        # Collect genotype matrices and frequencies
        for pop_name, pop in self.populations.items():
            results['populations'][pop_name] = {
                'genotype_matrix': pop.get_genotype_matrix(),
                'frequencies': pop.get_frequencies(),
                'history': pop.history
            }

        # Create truth table (which loci were selected and effect sizes)
        results['truth_table'] = {
            'selected_loci': selected_loci,
            'effect_sizes': {locus: 0.1 for locus in selected_loci},  # Fixed effect size for now
            'scenario_type': config['type']
        }

        return results

def save_results(results: Dict[str, Any], output_dir: str) -> None:
    """Save simulation results to files"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Save genotype matrices
    for pop_name, pop_data in results['populations'].items():
        genotype_file = Path(output_dir) / f"{pop_name}_genotypes.csv"
        pd.DataFrame(pop_data['genotype_matrix']).to_csv(genotype_file, index=False)

    # Save truth table
    truth_file = Path(output_dir) / "truth_table.json"
    with open(truth_file, 'w') as f:
        json.dump(results['truth_table'], f, indent=2)

    # Save summary
    summary_file = Path(output_dir) / "summary.json"
    summary = {
        'scenario_type': results['scenario_type'],
        'seed': results['seed'],
        'n_populations': len(results['populations']),
        'selected_loci': results['selected_loci']
    }
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)

def load_config(config_file: str) -> Dict[str, Any]:
    """Load configuration from YAML file"""
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def main():
    """Main function to run simulations from command line"""
    parser = argparse.ArgumentParser(description='Demographic Simulator')
    parser.add_argument('config_file', help='Path to configuration file')
    parser.add_argument('output_dir', help='Output directory for results')
    parser.add_argument('--scenario', help='Specific scenario to run (default: all)')

    args = parser.parse_args()

    # Load configuration
    config = load_config(args.config_file)

    # Create simulator
    simulator = DemographicSimulator(config)
    simulator.setup_populations()

    # Run scenarios
    scenarios = config.get('scenarios', {})

    if args.scenario and args.scenario in scenarios:
        # Run specific scenario
        scenario_config = scenarios[args.scenario]
        scenario_config['seed'] = config.get('seed', 42)
        results = simulator.run_scenario(scenario_config)
        save_results(results, Path(args.output_dir) / args.scenario)
    else:
        # Run all scenarios
        for scenario_name, scenario_config in scenarios.items():
            scenario_config['seed'] = config.get('seed', 42)
            print(f"Running scenario: {scenario_name}")
            results = simulator.run_scenario(scenario_config)
            save_results(results, Path(args.output_dir) / scenario_name)

if __name__ == "__main__":
    main()