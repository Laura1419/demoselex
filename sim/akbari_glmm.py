#!/usr/bin/env python3
"""
Akbari-style GLMM statistic implementation

Implementation of the GLMM from Akbari et al. 2026: allele frequency ~ time with
population random effects / GRM, leave-one-chromosome-out, locus-specific variance.

Key features:
- Binomial likelihood with logit link
- Population structure captured via GRM (Genetic Relatedness Matrix)
- Leave-one-chromosome-out cross-validation
- Produces Z-scores or posterior probabilities of selection
"""

import numpy as np
import pandas as pd
from scipy import linalg, stats
from scipy.optimize import minimize
from typing import Optional, Tuple, Dict, List
import warnings


class AkbariGLMM:
    """Akbari-style GLMM for detecting selection in ancient DNA time-series data.

    Implements the GLMM from Akbari et al. 2026:
    - Allele frequency ~ time with population random effects
    - Genetic Relatedness Matrix (GRM) for population structure
    - Leave-one-chromosome-out cross-validation
    - Produces Z-scores measuring selection strength

    Parameters
    ----------
    time_col : str
        Name of the column containing time values (years BP)
    allele_count_col : str
        Name of the column containing alternative allele counts
    total_count_col : str
        Name of the column containing total allele counts
    random_state : int, optional
        Random seed for reproducibility
    """

    def __init__(self,
                 time_col: str = "time_bp",
                 allele_count_col: str = "alt_count",
                 total_count_col: str = "total_count",
                 random_state: Optional[int] = None):

        self.time_col = time_col
        self.allele_count_col = allele_count_col
        self.total_count_col = total_count_col
        self.random_state = random_state

        # Model parameters (to be fitted)
        self.beta_time = None  # Selection coefficient (main parameter of interest)
        self.beta_intercept = None
        self.tau1 = None  # Variance component for GRM
        self.tau2 = None  # Variance component for identity matrix
        self.converged = False

        # Cache for computational efficiency
        self.V_inv = None
        self.K = None  # Genetic Relatedness Matrix
        self.X = None  # Design matrix

    def _check_inputs(self, data: pd.DataFrame):
        """Validate input data structure."""
        required_cols = [self.time_col, self.allele_count_col, self.total_count_col]
        missing_cols = [col for col in required_cols if col not in data.columns]

        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        if np.any(data[self.total_count_col] <= 0):
            raise ValueError("Total allele counts must be positive")

        if np.any(data[self.allele_count_col] < 0):
            raise ValueError("Allele counts cannot be negative")

        if np.any(data[self.allele_count_col] > data[self.total_count_col]):
            raise ValueError("Allele counts cannot exceed total counts")

    def compute_grm(self, genotypes: np.ndarray, method: str = "vanraden") -> np.ndarray:
        """Compute Genetic Relatedness Matrix from genotype data.

        Parameters
        ----------
        genotypes : np.ndarray
            Genotype matrix of shape (n_individuals, n_snps) with values in [0,1,2]
        method : str
            Method for GRM computation, either 'vanraden' or 'standardized'

        Returns
        -------
        K : np.ndarray
            Genetic Relatedness Matrix of shape (n_individuals, n_individuals)
        """
        n_individuals, n_snps = genotypes.shape

        if method == "vanraden":
            # VanRaden method (standard in genomic prediction)
            # Center genotypes by their mean
            P = genotypes.mean(axis=1, keepdims=True) / 2
            Z = genotypes - 2 * P

            # Compute GRM: Z Z^T / sum(2p(1-p))
            scaling = np.sum(2 * P * (1 - P), axis=1)
            K = Z @ Z.T / scaling.sum()

        elif method == "standardized":
            # Standardized GRM
            Z = (genotypes - genotypes.mean(axis=0)) / genotypes.std(axis=0, ddof=1)
            K = Z @ Z.T / n_snps

        else:
            raise ValueError("Method must be 'vanraden' or 'standardized'")

        # Ensure symmetry and positive semi-definiteness
        K = (K + K.T) / 2

        # Add small regularization for numerical stability
        epsilon = 1e-6
        K += epsilon * np.eye(n_individuals)

        return K

    def _design_matrix(self, data: pd.DataFrame) -> np.ndarray:
        """Construct design matrix with intercept and time columns."""
        n_samples = len(data)
        X = np.ones((n_samples, 2))  # intercept + time
        X[:, 1] = data[self.time_col].values
        return X

    def _working_response(self, y: np.ndarray, n: np.ndarray, p_hat: np.ndarray) -> np.ndarray:
        """Compute working response for PQL iteration."""
        # For binomial model with logit link
        return np.log(p_hat / (1 - p_hat)) + (y - n * p_hat) / (n * p_hat * (1 - p_hat) + 1e-8)

    def _variance_components_ai(self, X: np.ndarray, y_working: np.ndarray,
                               random_effects: np.ndarray, K: np.ndarray) -> Tuple[float, float]:
        """Estimate variance components using Average Information algorithm."""

        n = len(random_effects)

        # Current variance matrix
        V = self.tau1 * K + self.tau2 * np.eye(n)
        V_inv = linalg.inv(V)

        # Projection matrix
        XViX = X.T @ V_inv @ X
        XViX_inv = linalg.inv(XViX)
        P = V_inv - V_inv @ X @ XViX_inv @ X.T @ V_inv

        # Average Information matrix
        PK = P @ K
        PI = P @ np.eye(n)

        AI = np.array([
            [np.trace(PK @ PK), np.trace(PK @ PI)],
            [np.trace(PI @ PK), np.trace(PI @ PI)]
        ])

        # Score vector
        Py = P @ y_working
        score = np.array([
            random_effects.T @ K @ Py - np.trace(PK),
            random_effects.T @ Py - np.trace(PI)
        ])

        # Update variance components
        delta = linalg.solve(AI, score)
        tau1_new = max(1e-6, self.tau1 + delta[0])
        tau2_new = max(1e-6, self.tau2 + delta[1])

        return tau1_new, tau2_new

    def fit(self, data: pd.DataFrame, K: Optional[np.ndarray] = None,
            genotypes: Optional[np.ndarray] = None, max_iter: int = 100,
            tol: float = 1e-6) -> Dict:
        """Fit the GLMM model to data.

        Parameters
        ----------
        data : pd.DataFrame
            DataFrame containing time, allele counts, and total counts
        K : np.ndarray, optional
            Precomputed GRM matrix
        genotypes : np.ndarray, optional
            Genotype matrix for computing GRM if K not provided
        max_iter : int
            Maximum number of PQL iterations
        tol : float
            Convergence tolerance

        Returns
        -------
        dict
            Dictionary containing fitted parameters and convergence info
        """
        self._check_inputs(data)

        n_samples = len(data)
        y = data[self.allele_count_col].values
        n_total = data[self.total_count_col].values

        # Compute or use provided GRM
        if K is None:
            if genotypes is None:
                raise ValueError("Either K or genotypes must be provided")
            self.K = self.compute_grm(genotypes)
        else:
            self.K = K

        if self.K.shape != (n_samples, n_samples):
            raise ValueError("GRM shape must match number of samples")

        # Design matrix
        self.X = self._design_matrix(data)

        # Initialize parameters
        self.tau1 = 0.1  # initial GRM variance
        self.tau2 = 0.1  # initial residual variance

        # Initial fixed effects (logistic regression ignoring structure)
        p0 = np.clip(y / n_total, 0.01, 0.99)
        y_working0 = self._working_response(y, n_total, p0)

        # Initial GLS estimate (ignoring binomial nature for simplicity)
        V0 = self.tau1 * self.K + self.tau2 * np.eye(n_samples)
        V0_inv = linalg.inv(V0)
        XViX0 = self.X.T @ V0_inv @ self.X
        self.beta = linalg.solve(XViX0, self.X.T @ V0_inv @ y_working0)

        # PQL iterations
        prev_log_likelihood = -np.inf

        for iteration in range(max_iter):
            # Current linear predictor
            eta = self.X @ self.beta
            p_hat = 1 / (1 + np.exp(-eta))

            # Working response
            y_working = self._working_response(y, n_total, p_hat)

            # Update variance components using Average Information
            random_effects = y_working - eta
            tau1_new, tau2_new = self._variance_components_ai(
                self.X, y_working, random_effects, self.K
            )

            # Update fixed effects
            V_new = tau1_new * self.K + tau2_new * np.eye(n_samples)
            V_inv_new = linalg.inv(V_new)
            XViX_new = self.X.T @ V_inv_new @ self.X
            beta_new = linalg.solve(XViX_new, self.X.T @ V_inv_new @ y_working)

            # Check convergence
            beta_change = np.max(np.abs(beta_new - self.beta))
            tau_change = max(abs(tau1_new - self.tau1), abs(tau2_new - self.tau2))

            # Update parameters
            self.beta = beta_new
            self.tau1 = tau1_new
            self.tau2 = tau2_new
            self.V_inv = V_inv_new

            if max(beta_change, tau_change) < tol:
                self.converged = True
                break

        # Extract parameters
        self.beta_intercept = self.beta[0]
        self.beta_time = self.beta[1]

        # Compute standard errors and Z-score for beta_time
        XViX = self.X.T @ self.V_inv @ self.X
        beta_cov = linalg.inv(XViX)
        beta_time_se = np.sqrt(beta_cov[1, 1])
        z_score = self.beta_time / beta_time_se

        results = {
            'beta_time': self.beta_time,
            'beta_intercept': self.beta_intercept,
            'tau1': self.tau1,
            'tau2': self.tau2,
            'z_score': z_score,
            'se_beta_time': beta_time_se,
            'converged': self.converged,
            'n_iterations': iteration + 1
        }

        return results

    def leave_one_chromosome_out(self, data: pd.DataFrame, chromosome_col: str,
                                genotypes: np.ndarray, chromosome_ids: List[str]) -> Dict[str, Dict]:
        """Perform leave-one-chromosome-out cross-validation.

        Parameters
        ----------
        data : pd.DataFrame
            Full dataset
        chromosome_col : str
            Column name containing chromosome identifiers
        genotypes : np.ndarray
            Full genotype matrix
        chromosome_ids : List[str]
            List of chromosome IDs to leave out

        Returns
        -------
        Dict[str, Dict]
            Dictionary mapping chromosome IDs to their GLMM results
        """
        results = {}

        for chrom in chromosome_ids:
            # Split data into train (all but this chromosome) and test (this chromosome)
            train_mask = data[chromosome_col] != chrom
            test_mask = data[chromosome_col] == chrom

            if not np.any(test_mask):
                warnings.warn(f"No data for chromosome {chrom}, skipping")
                continue

            # Fit model on training data
            train_data = data[train_mask]
            train_genotypes = genotypes[train_mask]

            train_results = self.fit(train_data, genotypes=train_genotypes)

            if not train_results['converged']:
                warnings.warn(f"Model did not converge for chromosome {chrom} training")
                continue

            # Apply to test chromosome
            test_data = data[test_mask]

            # For test chromosome, we compute expected allele frequency
            # using the training model parameters
            X_test = self._design_matrix(test_data)
            eta_test = X_test @ np.array([train_results['beta_intercept'],
                                         train_results['beta_time']])
            p_test = 1 / (1 + np.exp(-eta_test))

            # Compute residual for test data
            y_test = test_data[self.allele_count_col].values
            n_test = test_data[self.total_count_col].values
            residuals = (y_test - n_test * p_test) / np.sqrt(n_test * p_test * (1 - p_test) + 1e-8)

            results[chrom] = {
                'train_results': train_results,
                'test_mean_residual': np.mean(residuals),
                'test_std_residual': np.std(residuals),
                'test_z_score': np.mean(residuals) / (np.std(residuals) + 1e-8)
            }

        return results


def simulate_data(n_individuals: int = 100, n_timepoints: int = 10,
                 selection_strength: float = 0.0, population_structure: bool = True,
                 random_state: Optional[int] = None) -> Tuple[pd.DataFrame, np.ndarray]:
    """Simulate data for testing the GLMM.

    Parameters
    ----------
    n_individuals : int
        Number of individuals
    n_timepoints : int
        Number of distinct time points
    selection_strength : float
        Strength of selection (beta_time parameter)
    population_structure : bool
        Whether to include population structure in simulation
    random_state : int, optional
        Random seed

    Returns
    -------
    Tuple[pd.DataFrame, np.ndarray]
        Simulated data and genotype matrix
    """
    rng = np.random.default_rng(random_state)

    # Simulate time points (spread over 10,000 years)
    times = np.linspace(0, 10000, n_timepoints)

    # Assign individuals to time points
    time_assignments = rng.choice(times, size=n_individuals)

    # Simulate genotypes with population structure
    if population_structure:
        # Create population clusters
        n_clusters = 3
        cluster_centers = rng.normal(0, 1, (n_clusters, 100))  # 100 SNP basis

        # Assign individuals to clusters
        cluster_probs = rng.dirichlet([1, 1, 1], size=n_individuals)
        cluster_assignments = np.argmax(cluster_probs, axis=1)

        # Simulate genotypes with cluster structure
        genotypes = np.zeros((n_individuals, 1000))  # 1000 SNPs
        for i in range(n_individuals):
            center = cluster_centers[cluster_assignments[i]]
            genotypes[i, :100] = rng.normal(center, 0.5, 100)
            genotypes[i, 100:] = rng.normal(0, 1, 900)

        # Binarize genotypes (convert to 0,1,2)
        genotypes = np.round(np.clip(genotypes, -3, 3)).astype(int)
        genotypes = np.where(genotypes < 0, 0, genotypes)
        genotypes = np.where(genotypes > 2, 2, genotypes)
    else:
        # No population structure
        genotypes = rng.binomial(2, 0.3, (n_individuals, 1000))

    # Simulate allele frequencies with selection
    intercept = -1.0  # Baseline log-odds

    # Linear predictor including selection
    eta = intercept + selection_strength * (time_assignments / 10000)
    p_true = 1 / (1 + np.exp(-eta))

    # Sample allele counts
    total_counts = rng.poisson(50, n_individuals)  # Sequencing depth
    alt_counts = rng.binomial(total_counts, p_true)

    # Create dataframe
    data = pd.DataFrame({
        'time_bp': time_assignments,
        'alt_count': alt_counts,
        'total_count': total_counts,
        'individual_id': range(n_individuals)
    })

    return data, genotypes


def test_calibration():
    """Test that the GLMM is calibrated under no selection (beta_time = 0)."""
    print("Testing calibration under no selection...")

    # Simulate 100 datasets under no selection
    n_datasets = 100
    z_scores = []

    for i in range(n_datasets):
        data, genotypes = simulate_data(selection_strength=0.0, random_state=i)

        glmm = AkbariGLMM()
        results = glmm.fit(data, genotypes=genotypes)

        if results['converged']:
            z_scores.append(results['z_score'])

    z_scores = np.array(z_scores)

    # Test calibration: z-scores should follow N(0,1) under null
    mean_z = np.mean(z_scores)
    std_z = np.std(z_scores)

    print(f"Mean Z-score: {mean_z:.3f} (expected: 0.0)")
    print(f"Std Z-score: {std_z:.3f} (expected: 1.0)")

    # KS test for normality
    ks_stat, p_value = stats.kstest(z_scores, 'norm')
    print(f"KS test p-value: {p_value:.3f}")

    return len(z_scores) > 0 and abs(mean_z) < 0.2 and 0.8 < std_z < 1.2


def test_power():
    """Test that the GLMM has power to detect strong selection."""
    print("Testing power under strong selection...")

    # Simulate 50 datasets with strong selection
    n_datasets = 50
    z_scores = []

    for i in range(n_datasets):
        data, genotypes = simulate_data(selection_strength=2.0, random_state=i + 1000)

        glmm = AkbariGLMM()
        results = glmm.fit(data, genotypes=genotypes)

        if results['converged']:
            z_scores.append(results['z_score'])

    z_scores = np.array(z_scores)

    # Under strong selection, z-scores should be significantly different from 0
    mean_z = np.mean(z_scores)
    significant = np.sum(np.abs(z_scores) > 1.96) / len(z_scores)

    print(f"Mean Z-score: {mean_z:.3f}")
    print(f"Proportion significant (|Z| > 1.96): {significant:.3f}")

    return len(z_scores) > 0 and significant > 0.8


if __name__ == "__main__":
    # Run tests
    print("Akbari-style GLMM Implementation Tests")
    print("=" * 50)

    calibration_ok = test_calibration()
    print(f"Calibration test: {'PASS' if calibration_ok else 'FAIL'}")
    print()

    power_ok = test_power()
    print(f"Power test: {'PASS' if power_ok else 'FAIL'}")
    print()

    # Example usage
    print("Example usage:")
    data, genotypes = simulate_data(n_individuals=200, selection_strength=1.0, random_state=42)

    glmm = AkbariGLMM()
    results = glmm.fit(data, genotypes=genotypes)

    print(f"Fitted beta_time (selection coefficient): {results['beta_time']:.4f}")
    print(f"Z-score: {results['z_score']:.4f}")
    print(f"Converged: {results['converged']}")

    if calibration_ok and power_ok:
        print("\nAll tests passed!")
    else:
        print("\nSome tests failed - check implementation")