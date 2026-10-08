"""Focused deterministic regressions for matched-locus rank calibration."""

import numpy as np
import pandas as pd
import pytest

from empirical_null import FEATURES, empirical_null
from trajectory_stat import trajectory_statistic


COLUMNS = ["locus", "trajectory", "time", "alt_count", "total_count"]


def _fixture():
    loci = ["candidate", "good_a", "good_b", "near", "other_chr", "other_time"]
    loci += [f"bad_{feature}" for feature in FEATURES]
    rows = []
    meta = []
    for i, locus in enumerate(loci):
        record = dict(locus=locus, chromosome="1", position=2_000_000 * i,
                      allele_frequency=0.3, recombination_rate=1.0, local_ld=0.2,
                      imputation_quality=0.95, missingness=0.03,
                      temporal_coverage=0.9)
        if locus == "near":
            record["position"] = 10
        if locus == "other_chr":
            record["chromosome"] = "2"
        if locus.startswith("bad_"):
            feature = locus[4:]
            record[feature] += 2 * {"allele_frequency": .05, "recombination_rate": .5,
                                   "local_ld": .1, "imputation_quality": -.05,
                                   "missingness": .05, "temporal_coverage": -.05}[feature]
        meta.append(record)
        for trajectory in ("a", "b", "c"):
            for time in ((1, 2, 3) if locus == "other_time" else (0, 1, 2)):
                # All-zero slopes are avoided, while score ties exercise the +1 rule.
                alt = 25 + 5 * time if locus != "good_b" else 30 - 4 * time
                rows.append((locus, trajectory, time, alt, 100))
    return pd.DataFrame(rows, columns=COLUMNS), pd.DataFrame(meta)


def test_all_matching_dimensions_and_support_are_enforced():
    counts, metadata = _fixture()
    result = empirical_null(counts, metadata, ["candidate"], min_controls=2)
    assert result.matches.control.tolist() == ["good_a", "good_b"]
    assert result.matches["rank"].tolist() == [1, 2]
    assert result.results.n_controls.iloc[0] == 2
    assert result.results.p_empirical.iloc[0] == pytest.approx(
        (1 + (result.matches.control_score >= result.results.replication_score.iloc[0]).sum()) / 3)
    assert {"concordance", "family_concordance", "z_score", "replication_score"} <= set(result.statistics)
    assert (result.matches.distance >= 0).all()
    # Input order never changes ties or matching; the outcome is not a matching feature.
    reordered = empirical_null(counts.iloc[::-1], metadata.iloc[::-1], ["candidate"], min_controls=2)
    pd.testing.assert_frame_equal(result.results, reordered.results)
    pd.testing.assert_frame_equal(result.matches, reordered.matches)


def test_nearest_controls_and_prespecified_family_covariance_are_reused():
    counts, metadata = _fixture()
    metadata.loc[metadata.locus == "good_b", "allele_frequency"] += .02
    rho = pd.DataFrame([[1, .3], [.3, 1]], index=["ab", "c"], columns=["ab", "c"])
    family = {"a": "ab", "b": "ab", "c": "c"}
    result = empirical_null(counts, metadata, ["candidate"], family=family,
                            relatedness=rho, min_controls=1, max_controls=1)
    assert result.matches.control.tolist() == ["good_a"]
    fitted = trajectory_statistic(counts, family=family, relatedness=rho).summary
    for locus in ("candidate", "good_a"):
        actual = result.statistics.set_index("locus").loc[locus]
        expected = fitted.set_index("locus").loc[locus]
        assert actual.replication_score == pytest.approx(
            abs(expected.z_score) * expected.family_concordance)
    with pytest.raises(ValueError, match="at least two prespecified"):
        empirical_null(counts, metadata, ["candidate"],
                       family={"a": "same", "b": "same", "c": "same"}, min_controls=1)


def test_conservative_ties_and_holdout_candidates():
    counts, metadata = _fixture()
    # A second candidate is removed from the background pool, not offered as
    # a matching control even if it has identical metadata and time support.
    both = empirical_null(counts, metadata, ["candidate", "good_a"], min_controls=1)
    assert "good_a" not in both.matches.query("locus == 'candidate'").control.to_list()
    assert "candidate" not in both.matches.query("locus == 'good_a'").control.to_list()
    tied = empirical_null(counts.query("locus in ['candidate', 'good_a']"),
                            metadata.query("locus in ['candidate', 'good_a']"),
                            ["candidate"], min_controls=1)
    assert tied.results.p_empirical.iloc[0] == 1
    assert tied.results.percentile.iloc[0] == 50


def test_fail_closed_on_no_controls_invalid_metadata_or_undefined_score():
    counts, metadata = _fixture()
    with pytest.raises(ValueError, match="only 2 matched controls"):
        empirical_null(counts, metadata, ["candidate"], min_controls=3)
    with pytest.raises(ValueError, match="match exactly"):
        empirical_null(counts.query("locus != 'good_b'"), metadata, ["candidate"], min_controls=1)
    with pytest.raises(ValueError, match="allele_frequency must be in"):
        empirical_null(counts, metadata.assign(allele_frequency=1.5), ["candidate"])
    with pytest.raises(ValueError, match="unknown matching caliper"):
        empirical_null(counts, metadata, ["candidate"], calipers={"unknown": 1})
    with pytest.raises(ValueError, match="defined, finite replication score"):
        constant = counts.copy()
        constant["alt_count"] = 20
        empirical_null(constant, metadata, ["candidate"], min_controls=1)
