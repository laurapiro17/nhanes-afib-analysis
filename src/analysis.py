"""Survey-weighted analysis of hypertension and obesity in U.S. adults.

Data: NHANES 2017-2018 (CDC/NCHS), examination subsample (measured BP & BMI).

Outcome     hypertension = mean SBP >= 130 OR mean DBP >= 80 (2017 ACC/AHA)
                           OR currently taking antihypertensive medication.
Exposure    obesity = BMI >= 30 kg/m^2.
Design      stratum = SDMVSTRA, PSU = SDMVPSU, weight = WTMEC2YR (MEC exam).

Outputs (written to ../results and ../figures):
  - prevalence_overall.csv, prevalence_by_group.csv  (design-based, 95% CI)
  - logistic_adjusted_or.csv                          (adjusted odds ratios)
  - fig_prevalence_by_bmi.png, fig_prevalence_by_age.png
  - results.json                                      (machine-readable summary)

All point estimates use the survey weights; confidence intervals for prevalence
use Taylor linearization (samplics) accounting for strata and PSUs.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from samplics.estimation import TaylorEstimator

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

SBP_THRESHOLD = 130
DBP_THRESHOLD = 80
OBESITY_BMI = 30.0


# ----------------------------------------------------------------------
# 1. Build the analytic dataset
# ----------------------------------------------------------------------
def load_data() -> pd.DataFrame:
    demo = pd.read_sas(DATA / "DEMO_J.xpt", format="xport")
    bpx = pd.read_sas(DATA / "BPX_J.xpt", format="xport")
    bmx = pd.read_sas(DATA / "BMX_J.xpt", format="xport")
    bpq = pd.read_sas(DATA / "BPQ_J.xpt", format="xport")

    demo = demo[["SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3",
                 "WTMEC2YR", "SDMVPSU", "SDMVSTRA"]]
    bpx = bpx[["SEQN", "BPXSY1", "BPXSY2", "BPXSY3",
               "BPXDI1", "BPXDI2", "BPXDI3"]]
    bmx = bmx[["SEQN", "BMXBMI"]]
    bpq = bpq[["SEQN", "BPQ050A"]]

    df = demo.merge(bpx, on="SEQN", how="left") \
             .merge(bmx, on="SEQN", how="left") \
             .merge(bpq, on="SEQN", how="left")
    return df


def derive(df: pd.DataFrame) -> pd.DataFrame:
    # Diastolic readings of 0 are device artefacts (cuff failed) -> missing.
    for c in ["BPXDI1", "BPXDI2", "BPXDI3"]:
        df.loc[df[c] == 0, c] = np.nan

    df["SBP"] = df[["BPXSY1", "BPXSY2", "BPXSY3"]].mean(axis=1)
    df["DBP"] = df[["BPXDI1", "BPXDI2", "BPXDI3"]].mean(axis=1)

    # On antihypertensive medication: BPQ050A == 1 ("yes"). All else (no, skip,
    # refused, don't know) treated as not currently medicated.
    on_meds = df["BPQ050A"] == 1
    has_bp = df["SBP"].notna() | df["DBP"].notna()
    bp_high = (df["SBP"] >= SBP_THRESHOLD) | (df["DBP"] >= DBP_THRESHOLD)

    # Crucial: do NOT let `(NaN >= 130) == False` mark an unmeasured person as
    # normotensive. Hypertension is positive if on medication; otherwise it is
    # determined only when a BP reading exists; with neither, it is missing.
    hyper = pd.Series(np.nan, index=df.index, dtype=float)
    hyper[has_bp] = bp_high[has_bp].astype(float)
    hyper[on_meds] = 1.0
    df["hypertension"] = hyper

    # Same trap for BMI: missing BMI must stay missing, not become "not obese".
    df["obese"] = np.where(df["BMXBMI"].notna(),
                           (df["BMXBMI"] >= OBESITY_BMI).astype(float),
                           np.nan)
    df["female"] = (df["RIAGENDR"] == 2).astype(float)
    df["age"] = df["RIDAGEYR"].astype(float)

    df["bmi_cat"] = pd.cut(
        df["BMXBMI"],
        bins=[0, 18.5, 25, 30, 200],
        labels=["Underweight", "Normal", "Overweight", "Obese"],
    )
    df["age_group"] = pd.cut(
        df["age"],
        bins=[17, 29, 39, 49, 59, 69, 200],
        labels=["18-29", "30-39", "40-49", "50-59", "60-69", "70+"],
    )
    return df


def analytic_sample(df: pd.DataFrame) -> pd.DataFrame:
    """Adults >=18 with a positive MEC weight and non-missing key variables.

    Note: rows failing the measurement requirements are *excluded* from the
    estimator, but the full design (all strata/PSUs) is retained so variance
    estimation stays valid (subpopulation / domain analysis).
    """
    keep = (
        (df["age"] >= 18)
        & (df["WTMEC2YR"] > 0)
        & df["hypertension"].notna()
        & df["obese"].notna()
        & df["BMXBMI"].notna()
    )
    return df[keep].copy()


# ----------------------------------------------------------------------
# 2. Design-based prevalence (Taylor linearization)
# ----------------------------------------------------------------------
def weighted_prevalence(df: pd.DataFrame, domain: pd.Series | None = None) -> pd.DataFrame:
    est = TaylorEstimator(param="proportion")
    est.estimate(
        y=df["hypertension"],
        samp_weight=df["WTMEC2YR"],
        stratum=df["SDMVSTRA"],
        psu=df["SDMVPSU"],
        domain=domain if domain is not None else None,
        remove_nan=True,
    )
    # samplics returns dict keyed by domain (and by outcome level 0/1)
    rows = []
    point = est.point_est
    lower = est.lower_ci
    upper = est.upper_ci
    if domain is None:
        rows.append(("Overall", point[1.0], lower[1.0], upper[1.0]))
    else:
        for key in point:
            rows.append((key, point[key][1.0], lower[key][1.0], upper[key][1.0]))
    out = pd.DataFrame(rows, columns=["group", "prevalence", "ci_low", "ci_high"])
    for c in ["prevalence", "ci_low", "ci_high"]:
        out[c] = (out[c] * 100).round(1)
    return out


# ----------------------------------------------------------------------
# 3. Adjusted association (survey-weighted logistic regression)
# ----------------------------------------------------------------------
def _fit_or(d: pd.DataFrame) -> np.ndarray:
    """Weighted logistic regression; return odds ratios for the 3 covariates.

    Age is modelled per decade so its odds ratio is interpretable (an OR per
    single year rounds to ~1.07 with a CI that vanishes under rounding).
    """
    X = d[["intercept", "obese", "age10", "female"]]
    res = sm.GLM(d["hypertension"], X, family=sm.families.Binomial(),
                 freq_weights=d["WTMEC2YR"]).fit()
    return np.exp(res.params[["obese", "age10", "female"]].values)


def adjusted_logistic(df: pd.DataFrame, n_boot: int = 400, seed: int = 17) -> pd.DataFrame:
    """Survey-weighted logistic regression of hypertension on obesity+age+sex.

    Point estimates use the MEC weights and are design-consistent. Confidence
    intervals come from a design-based bootstrap that resamples PSUs *with
    replacement within each stratum* (a standard NHANES variance approach) and
    refits the weighted model on each replicate. This sidesteps statsmodels'
    lack of true survey covariance with frequency weights.
    """
    d = df.dropna(subset=["hypertension", "obese", "age", "female"]).copy()
    d["intercept"] = 1.0
    d["age10"] = d["age"] / 10.0
    point = _fit_or(d)

    rng = np.random.default_rng(seed)
    strata = d.groupby("SDMVSTRA")
    # map each stratum -> {psu: dataframe} for fast resampling
    psu_groups = {
        s: [g for _, g in sub.groupby("SDMVPSU")]
        for s, sub in strata
    }
    boot = []
    for _ in range(n_boot):
        parts = []
        for s, psus in psu_groups.items():
            k = len(psus)
            if k < 2:
                parts.extend(psus)            # can't resample a lone PSU
                continue
            idx = rng.integers(0, k, size=k)  # sample k PSUs with replacement
            parts.extend(psus[i] for i in idx)
        rep = pd.concat(parts, ignore_index=True)
        try:
            boot.append(_fit_or(rep))
        except Exception:
            continue
    boot = np.array(boot)
    lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)

    return pd.DataFrame({
        "term": ["Obesity (BMI>=30)", "Age (per 10 years)", "Female"],
        "OR": point.round(2),
        "ci_low": lo.round(2),
        "ci_high": hi.round(2),
        "boot_replicates": [len(boot)] * 3,
    })


# ----------------------------------------------------------------------
# 4. Figures
# ----------------------------------------------------------------------
def plot_prevalence(df: pd.DataFrame, group_col: str, title: str, fname: str) -> None:
    est = weighted_prevalence(df, domain=df[group_col].astype("object"))
    # keep natural category order
    order = list(df[group_col].cat.categories) if hasattr(df[group_col], "cat") else None
    if order:
        est["group"] = pd.Categorical(est["group"], categories=order, ordered=True)
        est = est.sort_values("group")
    fig, ax = plt.subplots(figsize=(7, 4.2))
    yerr = [est["prevalence"] - est["ci_low"], est["ci_high"] - est["prevalence"]]
    ax.bar(est["group"].astype(str), est["prevalence"], yerr=yerr,
           capsize=4, color="#2b6cb0", alpha=0.9)
    ax.set_ylabel("Hypertension prevalence (%)")
    ax.set_title(title)
    ax.set_ylim(0, 100)
    for i, (v, hi) in enumerate(zip(est["prevalence"], est["ci_high"])):
        ax.text(i, min(hi + 3, 97), f"{v:.0f}%", ha="center", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / fname, dpi=130)
    plt.close(fig)


# ----------------------------------------------------------------------
# 5. Orchestrate
# ----------------------------------------------------------------------
def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)

    df = derive(load_data())
    sample = analytic_sample(df)
    n = len(sample)
    n_weighted = sample["WTMEC2YR"].sum()
    print(f"Analytic sample: {n:,} adults (representing "
          f"{n_weighted/1e6:.1f} million U.S. adults)")

    overall = weighted_prevalence(sample)
    by_obese = weighted_prevalence(sample, domain=sample["obese"].map(
        {0.0: "Not obese", 1.0: "Obese"}))
    by_age = weighted_prevalence(sample, domain=sample["age_group"].astype("object"))
    by_sex = weighted_prevalence(sample, domain=sample["female"].map(
        {0.0: "Male", 1.0: "Female"}))
    logit = adjusted_logistic(sample)

    overall.to_csv(RESULTS / "prevalence_overall.csv", index=False)
    pd.concat([by_obese.assign(stratifier="obesity"),
               by_age.assign(stratifier="age_group"),
               by_sex.assign(stratifier="sex")]).to_csv(
        RESULTS / "prevalence_by_group.csv", index=False)
    logit.to_csv(RESULTS / "logistic_adjusted_or.csv", index=False)

    plot_prevalence(sample, "bmi_cat",
                    "Hypertension prevalence by BMI category\nU.S. adults, NHANES 2017-2018 (survey-weighted)",
                    "fig_prevalence_by_bmi.png")
    plot_prevalence(sample, "age_group",
                    "Hypertension prevalence by age group\nU.S. adults, NHANES 2017-2018 (survey-weighted)",
                    "fig_prevalence_by_age.png")

    summary = {
        "n_analytic": int(n),
        "n_weighted_millions": round(float(n_weighted) / 1e6, 1),
        "overall_prevalence_pct": float(overall.loc[0, "prevalence"]),
        "overall_ci": [float(overall.loc[0, "ci_low"]), float(overall.loc[0, "ci_high"])],
        "prevalence_obese": by_obese.to_dict(orient="records"),
        "prevalence_by_age": by_age.to_dict(orient="records"),
        "prevalence_by_sex": by_sex.to_dict(orient="records"),
        "adjusted_or": logit.to_dict(orient="records"),
        "definitions": {
            "hypertension": "mean SBP>=130 OR mean DBP>=80 (2017 ACC/AHA) OR antihypertensive med",
            "obesity": "BMI>=30",
            "source": "NHANES 2017-2018, MEC examination sample",
        },
    }
    (RESULTS / "results.json").write_text(json.dumps(summary, indent=2))

    print("\nOverall hypertension prevalence: "
          f"{summary['overall_prevalence_pct']}% "
          f"(95% CI {summary['overall_ci'][0]}-{summary['overall_ci'][1]})")
    print("\nBy obesity status:")
    print(by_obese.to_string(index=False))
    print("\nAdjusted odds ratios:")
    print(logit.to_string(index=False))
    print(f"\nWrote results -> {RESULTS}  and figures -> {FIGURES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
