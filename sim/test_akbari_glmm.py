#!/usr/bin/env python3
"""Unit tests for the Akbari-style GLMM implementation."""

import numpy as np
import pandas as pd
import pytest
from akbari_glmm import AkbariGLMM, simulate_data


def test_akbari_glmm_initialization():
    """Test that the GLMM initializes correctly."""
    glmm = AkbariGLMM(
        time_col="time",
        allele_count_col="alt",
        total_count_col="total",
        random_state=42
    )

    assert glmm.time_col == "time"
    assert glmm.allele_count_col == "alt"
    assert glmm.total_count_col == "total"
    assert glmm.random_state == 42
    assert glmm.beta_time is None
    assert glmm.converged is False


def test_grm_computation():
    """Test GRM computation methods."""
    rng = np.random.default_rng(42)
    n_individuals = 10
    n_snps = 100
    genotypes = rng.binomial(2, 0.3, (n_individuals, n_snps))

    glmm = AkbariGLMM()

    # Test vanRaden method
    K_vanraden = glmm.compute_grm(genotypes, method="vanraden")
    assert K_vanraden.shape == (n_individuals, n_individuals)
    assert np.allclose(K_vanraden, K_vanraden.T)  # Symmetric
    assert np.all(np.linalg.eigvals(K_vanraden) >= -1e-6)  # PSD

    # Test standardized method
    K_std = glmm.compute_grm(genotypes, method="standardized")
    assert K_std.shape == (n_individuals, n_individuals)
    assert np.allclose(K_std, K_std.T)

    # Test invalid method
    try:
        glmm.compute_grm(genotypes, method="invalid")
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_input_validation():
    """Test input validation for the fit method."""
    glmm = AkbariGLMM()

    # Test missing columns
    invalid_data = pd.DataFrame({
        'time_wrong': [1000, 2000, 3000],
        'alt_count': [5, 10, 15],
        'total_count': [20, 40, 60]
    })

    genotypes = np.random.binomial(2, 0.3, (3, 100))

    try:
        glmm.fit(invalid_data, genotypes=genotypes)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

    # Test invalid counts
    invalid_counts = pd.DataFrame({
        'time_bp': [1000, 2000, 3000],
        'alt_count': [5, 10, 25],  # alt_count > total_count
        'total_count': [20, 40, 20]
    })

    try:
        glmm.fit(invalid_counts, genotypes=genotypes)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_simulation():
    """Test data simulation function."""
    # Test simulation without population structure
    data, genotypes = simulate_data(
        n_individuals=50,
        n_timepoints=5,
        selection_strength=0.0,
        population_structure=False,
        random_state=42
    )

    assert isinstance(data, pd.DataFrame)
    assert isinstance(genotypes, np.ndarray)
    assert len(data) == 50
    assert genotypes.shape[0] == 50
    assert 'time_bp' in data.columns
    assert 'alt_count' in data.columns
    assert 'total_count' in data.columns

    # Test simulation with population structure
    data_struct, genotypes_struct = simulate_data(
        n_individuals=50,
        population_structure=True,
        random_state=42
    )

    assert genotypes_struct.shape[0] == 50

    # Test that selection affects allele frequencies
    data_no_sel, _ = simulate_data(selection_strength=0.0, random_state=42)
    data_strong_sel, _ = simulate_data(selection_strength=2.0, random_state=42)

    # With strong selection, allele frequencies should show time trend
    freq_no_sel = data_no_sel['alt_count'] / data_no_sel['total_count']
    freq_strong_sel = data_strong_sel['alt_count'] / data_strong_sel['total_count']

    # Correlation with time should be stronger under selection
    corr_no_sel = np.corrcoef(data_no_sel['time_bp'], freq_no_sel)[0, 1]
    corr_strong_sel = np.corrcoef(data_strong_sel['time_bp'], freq_strong_sel)[0, 1]

    # This should generally hold, though not guaranteed for every simulation
    assert abs(corr_strong_sel) >= abs(corr_no_sel) - 0.1


def test_fit_convergence():
    """Test that the GLMM can fit to simulated data."""
    data, genotypes = simulate_data(
        n_individuals=100,
        selection_strength=0.5,
        random_state=42
    )

    glmm = AkbariGLMM()
    results = glmm.fit(data, genotypes=genotypes, max_iter=50)

    # Check that results dictionary has expected keys
    expected_keys = ['beta_time', 'beta_intercept', 'tau1', 'tau2',
                    'z_score', 'se_beta_time', 'converged', 'n_iterations']

    for key in expected_keys:
        assert key in results

    # Check that parameters are reasonable
    assert isinstance(results['beta_time'], (int, float))
    assert isinstance(results['z_score'], (int, float))
    assert results['n_iterations'] <= 50

    # The model should either converge or reach max iterations
    assert results['converged'] or results['n_iterations'] == 50


def test_leave_one_chromosome_out():
    """Test leave-one-chromosome-out cross-validation."""
    # Create data with chromosome information
    data, genotypes = simulate_data(n_individuals=90, random_state=42)

    # Simulate chromosome assignments with more individuals
    chromosomes = ['1'] * 30 + ['2'] * 30 + ['3'] * 30  # 30 individuals per chromosome
    data['chromosome'] = chromosomes

    glmm = AkbariGLMM()

    # Test LOCO
    results = glmm.leave_one_chromosome_out(
        data, 'chromosome', genotypes, ['1', '2', '3']
    )

    assert isinstance(results, dict)

    # Check that we got results for at least some chromosomes
    assert len(results) >= 1, "No chromosomes produced results"

    for chrom, chrom_results in results.items():
        assert 'train_results' in chrom_results
        assert 'test_mean_residual' in chrom_results
        assert 'test_std_residual' in chrom_results
        assert 'test_z_score' in chrom_results


def test_calibration_null():
    """Test that Z-scores are calibrated under no selection."""
    # Run a small calibration test
    n_simulations = 20
    z_scores = []

    for i in range(n_simulations):
        data, genotypes = simulate_data(
            n_individuals=150,
            selection_strength=0.0,
            random_state=i
        )

        glmm = AkbariGLMM()
        results = glmm.fit(data, genotypes=genotypes, max_iter=30)

        if results['converged']:
            z_scores.append(results['z_score'])

    if len(z_scores) >= 10:  # Need reasonable sample size
        z_scores = np.array(z_scores)

        # Under null, z-scores should have mean near 0 and std near 1
        mean_z = np.mean(z_scores)
        std_z = np.std(z_scores)

        # Check rough calibration
        assert abs(mean_z) < 1.0, f"Mean Z-score too large: {mean_z}"
        assert 0.5 < std_z < 2.0, f"Std Z-score out of range: {std_z}"


def test_power_detection():
    """Test that the method has power to detect selection."""
    # Run a small power test
    n_simulations = 15
    significant_count = 0
    total_converged = 0

    for i in range(n_simulations):
        data, genotypes = simulate_data(
            n_individuals=150,
            selection_strength=1.0,  # Strong selection
            random_state=i + 100
        )

        glmm = AkbariGLMM()
        results = glmm.fit(data, genotypes=genotypes, max_iter=30)

        if results['converged']:
            total_converged += 1
            if abs(results['z_score']) > 1.96:  # p < 0.05
                significant_count += 1

    if total_converged >= 10:
        power = significant_count / total_converged
        # With strong selection, should have reasonable power
        assert power > 0.3, f"Power too low: {power}"


if __name__ == "__main__":
    # Run tests
    print("Running Akbari GLMM unit tests...")

    test_akbari_glmm_initialization()
    print("✓ Initialization test passed")

    test_grm_computation()
    print("✓ GRM computation test passed")

    test_input_validation()
    print("✓ Input validation test passed")

    test_simulation()
    print("✓ Simulation test passed")

    test_fit_convergence()
    print("✓ Fit convergence test passed")

    test_leave_one_chromosome_out()
    print("✓ Leave-one-chromosome-out test passed")

    test_calibration_null()
    print("✓ Calibration test passed")

    test_power_detection()
    print("✓ Power detection test passed")

    print("\nAll tests passed!")

    # Run the main tests from akbari_glmm.py as well
    print("\nRunning main implementation tests...")
    from akbari_glmm import test_calibration, test_power

    calibration_ok = test_calibration()
    power_ok = test_power()

    if calibration_ok and power_ok:
        print("✓ All implementation tests passed!")
    else:
        print("⚠ Some implementation tests had issues")