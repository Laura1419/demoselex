"""Deterministic statistical and interface checks for the Akbari-style PQL model."""

import numpy as np
import pandas as pd
import pytest
from scipy.special import expit

from akbari_glmm import AkbariGLMM, simulate_data


def test_grm_centers_markers_not_people_and_loco_excludes_columns():
    rng = np.random.default_rng(42)
    g = rng.binomial(2, [0.2, 0.4, 0.5, 0.7], size=(12, 4)).astype(float)
    g[0, 0] = np.nan
    chrom = np.array(["1", "1", "2", "2"])
    k = AkbariGLMM.compute_grm(g, marker_chromosomes=chrom, target_chromosome="1")
    assert k.shape == (12, 12)  # chromosome does not index people
    assert np.allclose(k, k.T)
    assert np.linalg.eigvalsh(k).min() > -1e-10
    assert np.allclose(k.sum(axis=0), 0, atol=1e-10)
    changed = g.copy()
    changed[:, :2] = 2 - changed[:, :2]
    assert np.allclose(k, AkbariGLMM.compute_grm(
        changed, marker_chromosomes=chrom, target_chromosome="1"))
    assert not np.allclose(k, AkbariGLMM.compute_grm(
        changed, marker_chromosomes=chrom, target_chromosome="2"))
    assert np.isfinite(AkbariGLMM.compute_grm(g, method="standardized")).all()


def test_loco_loci_share_samples_and_exclude_target_chromosome():
    data, g = simulate_data(n_individuals=40, random_state=9, n_markers=60)
    chrom = np.repeat(["1", "2", "3"], 20)
    loci = {"neutral_1": data, "neutral_2": data.copy()}
    results = AkbariGLMM().leave_one_chromosome_out(
        loci, {"neutral_1": "1", "neutral_2": "2"}, g, chrom)
    for name, target in (("neutral_1", "1"), ("neutral_2", "2")):
        k = AkbariGLMM.compute_grm(g, marker_chromosomes=chrom,
                                    target_chromosome=target)
        direct = AkbariGLMM().fit(data, K=k)
        assert results[name]["z_score"] == pytest.approx(direct["z_score"])


def test_glm_limit_matches_binomial_fisher_information():
    data, g = simulate_data(n_individuals=56, random_state=40)
    r = AkbariGLMM().fit(data, genotypes=g, fixed_tau=(0, 0))
    t = data.time_bp.to_numpy()
    x = np.column_stack((np.ones(len(t)), t))
    p = expit(r["beta_intercept"] + r["beta_time"] * t)
    assert np.max(np.abs(x.T @ (data.alt_count - data.total_count * p))) < 1e-4
    info = x.T @ (data.total_count.to_numpy()[:, None] * p[:, None] *
                  (1 - p)[:, None] * x)
    assert r["se_beta_time"] == pytest.approx(np.sqrt(np.linalg.inv(info)[1, 1]), rel=1e-5)
    assert r["p_value"] < 1 and r["converged"]


def test_null_calibration_without_skipping_nonconvergence():
    # Data are generated from the fitted model with a structured u and 30
    # independent fixed seeds. This bounds gross inflation/bias, not genome-wide
    # tail calibration (which requires the downstream sim-sweeps study).
    scores = []
    for seed in range(30):
        d, g = simulate_data(n_individuals=56, random_state=seed)
        scores.append(AkbariGLMM().fit(d, genotypes=g)["z_score"])
    scores = np.asarray(scores)
    assert abs(scores.mean()) < 0.5
    assert 0.65 < scores.std(ddof=1) < 1.35
    assert np.mean(np.abs(scores) > 1.96) <= 0.15


def test_strong_single_locus_sweep_is_significant():
    # The background markers remain neutral; only tested counts have time trend.
    for seed in range(100, 110):
        data, g = simulate_data(n_individuals=56, selection_strength=1.2,
                                random_state=seed)
        result = AkbariGLMM().fit(data, genotypes=g)
        assert result["z_score"] > 5.45
        assert result["p_value"] < 5e-8


@pytest.mark.parametrize("strength,reference_z", [
    (0.0, -0.7681044214992059), (1.2, 8.142205576678593),
])
def test_released_pqlseqpy_fixed_tau_fixture(strength, reference_z):
    # Independent published PQLseqPy source commit 0934722, GLMM(X,Y,K,
    # fixed_tau=(.12,.04)).fit(), on this exact seeded fixture. X uses
    # standardized time and Y=[alt, total-alt]. Comparison was run separately;
    # the test has no runtime dependency on PQLseqPy/statsmodels/network.
    data, g = simulate_data(n_individuals=48, selection_strength=strength,
                            random_state=42)
    r = AkbariGLMM().fit(data, genotypes=g, fixed_tau=(0.12, 0.04))
    assert r["z_score"] == pytest.approx(reference_z, abs=0.003)


def test_invalid_inputs_fail_closed():
    d, g = simulate_data(n_individuals=20, random_state=8)
    with pytest.raises(ValueError, match="LOCO requires"):
        AkbariGLMM().fit(d, genotypes=g, target_chromosome="1")
    with pytest.raises(ValueError, match="exactly one"):
        AkbariGLMM().fit(d, K=np.eye(20), genotypes=g)
    with pytest.raises(ValueError, match="counts must"):
        AkbariGLMM().fit(d.assign(alt_count=-1), genotypes=g)
    with pytest.raises(ValueError, match="positive semidefinite"):
        AkbariGLMM().fit(d, K=-np.eye(20))
    with pytest.raises(ValueError, match="no background markers"):
        AkbariGLMM.compute_grm(g, marker_chromosomes=["1"] * g.shape[1],
                                target_chromosome="1")
    with pytest.raises(ValueError, match="seed"):
        simulate_data()
