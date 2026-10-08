"""Deterministic regression tests for trajectory-level replication statistics."""

import numpy as np
import pandas as pd
import pytest

from trajectory_stat import SLOPE_COLUMNS, trajectory_statistic


def _designed(seed=31415):
    rng = np.random.default_rng(seed)
    rows = []
    for locus, trend in (("up", 0.45), ("down", -0.45)):
        for trajectory, intercept in (("a", 0.3), ("b", 0.45), ("c", 0.6)):
            for time in np.linspace(0, 1, 9):
                depth = 2000
                probability = intercept + trend * (time - 0.5)
                rows.append((locus, trajectory, time,
                             rng.binomial(depth, probability), depth))
    return pd.DataFrame(rows, columns=("locus", "trajectory", "time",
                                       "alt_count", "total_count"))


def test_designed_signal_recovers_both_signs_and_tidy_schema():
    data = _designed()
    result = trajectory_statistic(data)
    assert set(SLOPE_COLUMNS).issubset(result.slopes.columns)
    assert len(result.slopes) == 6
    assert (result.slopes.n == 9).all()
    assert (result.slopes.span == 1).all()
    assert (result.slopes.se > 0).all()
    for locus, sign in (("up", 1), ("down", -1)):
        slope = result.slopes.query("locus == @locus")
        summary = result.summary.set_index("locus").loc[locus]
        assert (slope.direction == sign).all()
        assert (slope.slope * sign > 0.3).all()  # trend >> binomial noise
        assert summary.meta_slope * sign > 0.3
        assert summary.concordance == summary.family_concordance == 1
        assert summary.n_effective == 3
        assert summary.q_df == 2
    # Input row order, not a hidden random seed, determines no results.
    shuffled = data.sample(frac=1, random_state=77)
    pd.testing.assert_frame_equal(result.slopes, trajectory_statistic(shuffled).slopes)
    pd.testing.assert_frame_equal(result.summary, trajectory_statistic(shuffled).summary)


def test_weighted_slope_and_binomial_se_analytic_two_time_fixture():
    data = pd.DataFrame({"locus": ["x", "x"], "trajectory": ["a", "a"],
                         "time": [0, 2], "alt_count": [2, 8],
                         "total_count": [10, 10]})
    result = trajectory_statistic(data)
    row = result.slopes.iloc[0]
    assert row.slope == pytest.approx(0.3)
    p0, p1 = 2.5 / 11, 8.5 / 11
    assert row.se == pytest.approx(np.sqrt((p0 * (1 - p0) + p1 * (1 - p1)) / 10) / 2)
    assert row.n == 2 and row.span == 2 and row.direction == 1
    assert np.isnan(result.summary.q.iloc[0])
    assert np.isnan(result.summary.q_p_value.iloc[0])
    # Duplicate-time rows are combined as binomial trials, not time replicates.
    split = pd.DataFrame({"locus": ["x"] * 4, "trajectory": ["a"] * 4,
                          "time": [0, 0, 2, 2], "alt_count": [1, 1, 3, 5],
                          "total_count": [5] * 4})
    pd.testing.assert_frame_equal(result.slopes, trajectory_statistic(split).slopes)


def test_relatedness_and_family_collapse_discount_copies():
    data = _designed().query("locus == 'up'")
    independent = trajectory_statistic(data)
    labels = ["a", "b", "c"]
    rho = pd.DataFrame(np.full((3, 3), 0.4) + 0.6 * np.eye(3),
                       index=labels, columns=labels)
    correlated = trajectory_statistic(data, relatedness=rho)
    assert correlated.summary.n_effective.iloc[0] == pytest.approx(9 / 5.4)
    assert correlated.summary.meta_se.iloc[0] > independent.summary.meta_se.iloc[0]
    assert 0 <= correlated.summary.q_p_value.iloc[0] <= 1
    assert 0 <= correlated.summary.i_squared.iloc[0] <= 1 or np.isnan(correlated.summary.i_squared.iloc[0])

    copy = data.query("trajectory == 'a'").assign(trajectory="a_copy")
    duplicated = trajectory_statistic(pd.concat([data, copy], ignore_index=True),
                                      family={"a": "a", "a_copy": "a", "b": "b", "c": "c"},
                                      relatedness=rho)
    # Pooled Model 0 legitimately reweights the copied path; family-level
    # Model 2 evidence and effective count must not gain precision.
    ignored = ["n_trajectories", "pooled_slope", "pooled_se"]
    pd.testing.assert_series_equal(duplicated.summary.iloc[0].drop(ignored),
                                   correlated.summary.iloc[0].drop(ignored))
    assert duplicated.summary.n_trajectories.iloc[0] == 4


def test_turnover_only_can_create_pooled_model0_slope():
    rows = []
    # A and B remain at 0.2 and 0.8 while their abundance swaps.
    for time, depths in ((0, (90, 10)), (1, (10, 90))):
        for trajectory, depth, frequency in zip(("a", "b"), depths, (0.2, 0.8)):
            rows.append(("replacement", trajectory, time, int(depth * frequency), depth))
    data = pd.DataFrame(rows, columns=("locus", "trajectory", "time",
                                       "alt_count", "total_count"))
    result = trajectory_statistic(data)
    assert (result.slopes.slope == 0).all()
    assert np.isnan(result.summary.concordance.iloc[0])
    assert result.summary.pooled_slope.iloc[0] == pytest.approx(0.48)
    assert result.summary.meta_slope.iloc[0] == 0


def test_ties_zero_slopes_and_heterogeneity():
    rows = []
    for trajectory, (start, end) in {"a": (2, 8), "b": (8, 2), "c": (5, 5)}.items():
        for time, count in ((0, start), (1, end)):
            rows.append(("x", trajectory, time, count, 10))
    data = pd.DataFrame(rows, columns=("locus", "trajectory", "time",
                                       "alt_count", "total_count"))
    result = trajectory_statistic(data)
    assert dict(zip(result.slopes.trajectory, result.slopes.direction)) == {
        "a": 1, "b": -1, "c": 0}
    assert result.summary.concordance.iloc[0] == 0.5
    assert result.summary.q.iloc[0] > 0
    assert np.isnan(trajectory_statistic(data.query("trajectory == 'c'"))
                    .summary.concordance.iloc[0])


def test_invalid_data_and_correlation_fail_closed():
    data = _designed().query("locus == 'up'")
    with pytest.raises(ValueError, match="counts must"):
        trajectory_statistic(data.assign(alt_count=-1))
    with pytest.raises(ValueError, match="distinct times"):
        trajectory_statistic(data.query("time == 0"))
    with pytest.raises(ValueError, match="family must map"):
        trajectory_statistic(data, family={"a": "a"})
    labels = ["a", "b", "c"]
    singular = pd.DataFrame(np.ones((3, 3)), index=labels, columns=labels)
    with pytest.raises(ValueError, match="positive definite"):
        trajectory_statistic(data, relatedness=singular)
    negative = pd.DataFrame(np.eye(3), index=labels, columns=labels)
    negative.loc["a", "b"] = negative.loc["b", "a"] = -0.1
    with pytest.raises(ValueError, match="entries in"):
        trajectory_statistic(data, relatedness=negative)
