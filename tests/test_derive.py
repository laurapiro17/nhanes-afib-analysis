"""Unit tests for the analytic-variable derivation.

These run on small synthetic frames (no NHANES download needed), so they are
fast and deterministic in CI. They pin down the missing-data handling that a
naive analysis gets wrong — the difference between an unmeasured person being
*missing* versus being silently counted as a negative.
"""
import numpy as np
import pandas as pd

from analysis import derive, analytic_sample, OBESITY_BMI


def _raw(rows):
    """Build a raw-column frame like the merged NHANES extract derive() expects."""
    cols = ["SEQN", "RIDAGEYR", "RIAGENDR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA",
            "BPXSY1", "BPXSY2", "BPXSY3", "BPXDI1", "BPXDI2", "BPXDI3",
            "BMXBMI", "BPQ050A"]
    return pd.DataFrame(rows, columns=cols)


def test_high_bp_is_hypertensive_when_not_medicated():
    df = derive(_raw([[1, 50, 1, 1.0, 1, 1, 140, 142, 138, 88, 90, 86, 27.0, 2]]))
    assert df.loc[0, "hypertension"] == 1.0


def test_normal_bp_is_not_hypertensive():
    df = derive(_raw([[1, 50, 1, 1.0, 1, 1, 110, 112, 108, 70, 72, 68, 27.0, 2]]))
    assert df.loc[0, "hypertension"] == 0.0


def test_unmeasured_bp_stays_missing_not_zero():
    """The core trap: no BP reading and not on meds -> missing, NOT normotensive."""
    df = derive(_raw([[1, 50, 1, 1.0, 1, 1,
                       np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 27.0, 2]]))
    assert np.isnan(df.loc[0, "hypertension"])


def test_on_medication_is_hypertensive_even_without_high_bp():
    df = derive(_raw([[1, 50, 1, 1.0, 1, 1, 110, 112, 108, 70, 72, 68, 27.0, 1]]))
    assert df.loc[0, "hypertension"] == 1.0


def test_on_medication_overrides_missing_bp():
    df = derive(_raw([[1, 50, 1, 1.0, 1, 1,
                       np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 27.0, 1]]))
    assert df.loc[0, "hypertension"] == 1.0


def test_zero_diastolic_is_treated_as_artefact():
    """DBP readings of 0 are cuff failures; they must not pull the mean down."""
    df = derive(_raw([[1, 50, 1, 1.0, 1, 1, 110, 112, 108, 0, 0, 0, 27.0, 2]]))
    assert np.isnan(df.loc[0, "DBP"])
    # SBP normal, DBP unmeasurable -> still classified on SBP alone, not hypertensive
    assert df.loc[0, "hypertension"] == 0.0


def test_obesity_threshold_and_missing_bmi():
    df = derive(_raw([
        [1, 50, 1, 1.0, 1, 1, 110, 110, 110, 70, 70, 70, OBESITY_BMI, 2],       # exactly 30 -> obese
        [2, 50, 1, 1.0, 1, 1, 110, 110, 110, 70, 70, 70, 24.0, 2],              # < 30 -> not obese
        [3, 50, 1, 1.0, 1, 1, 110, 110, 110, 70, 70, 70, np.nan, 2],            # missing BMI -> missing
    ]))
    assert df.loc[0, "obese"] == 1.0
    assert df.loc[1, "obese"] == 0.0
    assert np.isnan(df.loc[2, "obese"])


def test_analytic_sample_excludes_minors_zero_weight_and_missing():
    df = derive(_raw([
        [1, 50, 1, 1.0, 1, 1, 140, 140, 140, 90, 90, 90, 31.0, 2],             # keep
        [2, 15, 1, 1.0, 1, 1, 140, 140, 140, 90, 90, 90, 31.0, 2],             # minor -> drop
        [3, 50, 1, 0.0, 1, 1, 140, 140, 140, 90, 90, 90, 31.0, 2],             # zero weight -> drop
        [4, 50, 1, 1.0, 1, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 2],  # all missing -> drop
    ]))
    keep = analytic_sample(df)
    assert list(keep["SEQN"]) == [1]
