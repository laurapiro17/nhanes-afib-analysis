"""Tests for the survey-weighted prevalence estimator on small synthetic designs."""
import numpy as np
import pandas as pd

from analysis import weighted_prevalence


def _design(hyper_values):
    """A minimal valid survey design: 2 strata x 2 PSUs, equal weights."""
    n = len(hyper_values)
    rng = np.arange(n)
    return pd.DataFrame({
        "hypertension": np.asarray(hyper_values, dtype=float),
        "WTMEC2YR": np.ones(n) * 1000.0,
        "SDMVSTRA": 1 + (rng // (n // 2)),     # two strata
        "SDMVPSU": 1 + ((rng // (n // 4)) % 2),  # two PSUs within each stratum
    })


def test_overall_prevalence_columns_and_range():
    df = _design([1, 0, 1, 0, 1, 0, 1, 0])
    out = weighted_prevalence(df)
    assert list(out.columns) == ["group", "prevalence", "ci_low", "ci_high"]
    assert len(out) == 1
    assert 0.0 <= out.loc[0, "prevalence"] <= 100.0


def test_all_positive_gives_full_prevalence():
    df = _design([1, 1, 1, 1, 1, 1, 1, 1])
    out = weighted_prevalence(df)
    assert out.loc[0, "prevalence"] == 100.0


def test_equal_weights_recover_sample_proportion():
    # 6 of 8 positive, equal weights -> 75%
    df = _design([1, 1, 1, 0, 1, 1, 0, 1])
    out = weighted_prevalence(df)
    assert out.loc[0, "prevalence"] == 75.0
