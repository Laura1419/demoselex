#!/usr/bin/env python3
"""
Test the demographic simulator with a simple scenario
"""

import sys
import os
import numpy as np

# Add parent directories to path
sys.path.append(os.path.dirname(__file__))
from simulator import DemographicSimulator

def test_neutral_drift():
    """Test neutral drift scenario"""
    print("Testing neutral drift scenario...")

    config = {
        'seed': 12345,
        'n_loci': 100,
        'populations': {
            'pop1': {'size': 100, 'time_origin': 0},
            'pop2': {'size': 100, 'time_origin': 0}
        }
    }

    simulator = DemographicSimulator(config)
    simulator.setup_populations()

    # Get initial frequencies
    initial_freqs = {}
    for pop_name, pop in simulator.populations.items():
        initial_freqs[pop_name] = pop.get_frequencies()

    # Run neutral drift
    scenario_config = {
        'type': 'neutral_drift',
        'generations': 10,
        'migration_rate': 0.01,
        'seed': 12345
    }

    results = simulator.run_scenario(scenario_config)

    # Check that frequencies changed (due to drift)
    final_freqs = {}
    for pop_name, pop in simulator.populations.items():
        final_freqs[pop_name] = pop.get_frequencies()

    # Verify that frequencies are different (drift occurred)
    for pop_name in initial_freqs:
        initial = initial_freqs[pop_name]
        final = final_freqs[pop_name]
        diff = np.mean(np.abs(initial - final))

        if diff > 0:
            print(f"✓ {pop_name}: Drift occurred (mean frequency change: {diff:.4f})")
        else:
            print(f"✗ {pop_name}: No drift detected")

    return results

def test_selection():
    """Test selection scenario"""
    print("\nTesting selection scenario...")

    config = {
        'seed': 12345,
        'n_loci': 100,
        'populations': {
            'pop1': {'size': 100, 'time_origin': 0}
        }
    }

    simulator = DemographicSimulator(config)
    simulator.setup_populations()

    # Get initial frequencies at selected locus
    selected_locus = 0
    initial_freq = simulator.populations['pop1'].get_frequencies()[selected_locus]

    # Run selection
    scenario_config = {
        'type': 'one_population_selection',
        'generations': 50,
        'selection_strength': 0.2,
        'selected_population': 'pop1',
        'selected_loci': [selected_locus],
        'seed': 12345
    }

    results = simulator.run_scenario(scenario_config)

    # Get final frequency
    final_freq = simulator.populations['pop1'].get_frequencies()[selected_locus]

    # Check if frequency increased (positive selection)
    if final_freq > initial_freq:
        print(f"✓ Selection occurred (frequency: {initial_freq:.4f} → {final_freq:.4f})")
    else:
        print(f"✗ Selection may not have worked (frequency: {initial_freq:.4f} → {final_freq:.4f})")

    return results

def test_truth_table():
    """Test that truth table is correctly generated"""
    print("\nTesting truth table generation...")

    config = {
        'seed': 12345,
        'n_loci': 100,
        'populations': {
            'pop1': {'size': 100, 'time_origin': 0}
        }
    }

    simulator = DemographicSimulator(config)
    simulator.setup_populations()

    selected_loci = [0, 10, 20]
    scenario_config = {
        'type': 'one_population_selection',
        'generations': 10,
        'selection_strength': 0.1,
        'selected_population': 'pop1',
        'selected_loci': selected_loci,
        'seed': 12345
    }

    results = simulator.run_scenario(scenario_config)

    # Check truth table
    truth_table = results['truth_table']

    if 'selected_loci' in truth_table and 'effect_sizes' in truth_table:
        print("✓ Truth table contains required fields")

        if sorted(truth_table['selected_loci']) == sorted(selected_loci):
            print("✓ Correct selected loci recorded")
        else:
            print("✗ Selected loci mismatch")
    else:
        print("✗ Missing fields in truth table")

    return results

def main():
    """Run all tests"""
    print("=== Testing Demographic Simulator ===\n")

    # Set numpy random seed for reproducibility
    np.random.seed(12345)

    # Run tests
    test_neutral_drift()
    test_selection()
    test_truth_table()

    print("\n=== Test Complete ===")

if __name__ == "__main__":
    main()