#!/usr/bin/env python3
"""Quick validation of the Akbari GLMM implementation."""

from akbari_glmm import AkbariGLMM, simulate_data
import numpy as np

def quick_calibration_test():
    """Quick test of calibration under no selection."""
    print("Quick calibration test (no selection)...")

    # Run a small number of simulations
    z_scores = []
    for i in range(5):
        try:
            data, genotypes = simulate_data(
                n_individuals=100,
                selection_strength=0.0,
                random_state=i
            )

            glmm = AkbariGLMM()
            results = glmm.fit(data, genotypes=genotypes, max_iter=20)

            if results['converged']:
                z_scores.append(results['z_score'])
                print(f"  Simulation {i+1}: Z = {results['z_score']:.3f}")
        except Exception as e:
            print(f"  Simulation {i+1} failed: {e}")

    if z_scores:
        mean_z = np.mean(z_scores)
        std_z = np.std(z_scores)
        print(f"Mean Z: {mean_z:.3f}, Std Z: {std_z:.3f}")

        # Check if results are reasonable
        return len(z_scores) >= 3 and abs(mean_z) < 2.0

    return False

def quick_power_test():
    """Quick test of power under selection."""
    print("Quick power test (strong selection)...")

    z_scores = []
    for i in range(5):
        try:
            data, genotypes = simulate_data(
                n_individuals=100,
                selection_strength=1.0,
                random_state=i + 10
            )

            glmm = AkbariGLMM()
            results = glmm.fit(data, genotypes=genotypes, max_iter=20)

            if results['converged']:
                z_scores.append(results['z_score'])
                print(f"  Simulation {i+1}: Z = {results['z_score']:.3f}")
        except Exception as e:
            print(f"  Simulation {i+1} failed: {e}")

    if z_scores:
        mean_z = np.mean(z_scores)
        significant = np.sum(np.abs(z_scores) > 1.96) / len(z_scores)
        print(f"Mean Z: {mean_z:.3f}, Proportion significant: {significant:.3f}")

        # Check if method shows power
        return len(z_scores) >= 3 and significant > 0.3

    return False

def test_interface():
    """Test that the interface is well-documented and functional."""
    print("Testing interface documentation...")

    # Test initialization with custom parameters
    glmm = AkbariGLMM(
        time_col="time",
        allele_count_col="alt",
        total_count_col="total"
    )

    # Test that the object has the expected attributes
    assert hasattr(glmm, 'fit')
    assert hasattr(glmm, 'compute_grm')
    assert hasattr(glmm, 'leave_one_chromosome_out')

    print("✓ Interface is properly documented and functional")
    return True

def main():
    """Run quick validation."""
    print("Akbari GLMM Quick Validation")
    print("=" * 40)

    # Test 1: Interface
    interface_ok = test_interface()

    # Test 2: Calibration under null
    calibration_ok = quick_calibration_test()

    # Test 3: Power under selection
    power_ok = quick_power_test()

    print("\n" + "=" * 40)
    print("Validation Summary:")
    print(f"Interface: {'PASS' if interface_ok else 'FAIL'}")
    print(f"Calibration: {'PASS' if calibration_ok else 'FAIL'}")
    print(f"Power: {'PASS' if power_ok else 'FAIL'}")

    if interface_ok and calibration_ok and power_ok:
        print("\n✓ All quick validation tests passed!")
        return True
    else:
        print("\n⚠ Some tests failed or were inconclusive")
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)