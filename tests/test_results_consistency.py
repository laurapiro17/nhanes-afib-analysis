"""Guards the committed headline results in results/results.json.

Catches a regression or an accidental hand-edit that would make the README's
numbers internally inconsistent (CIs not bracketing the point estimate,
prevalences out of range, or the reported obesity gradient reversing).
"""
import json
from pathlib import Path

import pytest

RESULTS = Path(__file__).resolve().parent.parent / "results" / "results.json"


@pytest.fixture(scope="module")
def results():
    return json.loads(RESULTS.read_text())


def test_overall_prevalence_in_range_and_bracketed(results):
    p = results["overall_prevalence_pct"]
    lo, hi = results["overall_ci"]
    assert 0.0 <= lo <= p <= hi <= 100.0


def test_obesity_gradient_holds(results):
    by = {r["group"]: r for r in results["prevalence_obese"]}
    assert by["Obese"]["prevalence"] > by["Not obese"]["prevalence"]
    for r in results["prevalence_obese"]:
        assert 0.0 <= r["ci_low"] <= r["prevalence"] <= r["ci_high"] <= 100.0


def test_weighted_population_is_plausible(results):
    # NHANES 2017-2018 represents ~230-250M U.S. adults; guard against a units slip.
    assert 200.0 <= results["n_weighted_millions"] <= 270.0
    assert results["n_analytic"] > 3000


def test_adjusted_obesity_or_above_one_and_bracketed(results):
    obesity = next(r for r in results["adjusted_or"] if "Obesity" in r["term"])
    assert obesity["ci_low"] > 1.0
    assert obesity["ci_low"] <= obesity["OR"] <= obesity["ci_high"]
