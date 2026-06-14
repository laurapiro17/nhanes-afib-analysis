# nhanes-afib-analysis

**A reproducible, survey-weighted analysis of hypertension and obesity in U.S.
adults** using public NHANES 2017–2018 data.

This is a from-scratch epidemiology pipeline: it downloads the raw CDC data,
builds the analytic dataset with correct missing-data handling, and produces
**design-based** estimates that respect the NHANES complex survey design
(stratification, clustering, and sampling weights) — not a naïve unweighted
summary.

[![CI](https://github.com/laurapiro17/nhanes-afib-analysis/actions/workflows/ci.yml/badge.svg)](https://github.com/laurapiro17/nhanes-afib-analysis/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.9-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Data](https://img.shields.io/badge/data-NHANES%202017--2018-orange)

## Headline result

Among an estimated **236 million U.S. adults**, hypertension prevalence
(2017 ACC/AHA definition) was **47.9% (95% CI 44.9–51.0)**, rising monotonically
with BMI and age. Obesity carried an adjusted odds ratio of **2.68 (2.33–3.08)**
after controlling for age and sex.

| | Hypertension prevalence |
|---|---|
| Not obese | 39.3% (36.3–42.4) |
| **Obese (BMI ≥ 30)** | **59.7% (56.0–63.4)** |

See [`REPORT.md`](REPORT.md) for the full write-up, methods, and figures.

## Why the methods matter

NHANES is a *complex survey*, not a simple random sample. Getting the numbers
right requires three things this repo does and a naïve analysis skips:

1. **Sampling weights** (`WTMEC2YR`) so estimates represent the U.S. population,
   not the (deliberately non-representative) sample.
2. **Design-based variance** — confidence intervals from Taylor linearization
   over strata (`SDMVSTRA`) and PSUs (`SDMVPSU`) via `samplics`, plus a
   **PSU bootstrap** for the regression. Ignoring clustering understates SEs.
3. **Honest missing data** — a person with no BP reading is *missing*, not
   "normotensive". `(NaN >= 130)` silently evaluates to `False`; coding the
   outcome naïvely would misclassify everyone unmeasured. The pipeline guards
   against exactly this.

## Reproduce it

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python src/fetch_data.py      # downloads 4 public NHANES files into data/
python src/analysis.py        # writes results/ and figures/, prints the summary
```

Runtime: ~20s (the 400-replicate bootstrap dominates).

## Tests

```bash
pip install pytest && pytest -q
```

The suite runs on synthetic frames (no data download needed), so CI is fast and
deterministic. It pins down the missing-data handling — that an unmeasured person
stays *missing* rather than being silently coded normotensive — the weighted
prevalence estimator on a small known design, and the internal consistency of the
committed `results.json` (CIs bracket the point estimates, obesity gradient holds).

## What's in here

```
src/fetch_data.py   download the NHANES public files (cached)
src/analysis.py     dataset construction, design-based prevalence, bootstrap logistic regression, figures
results/            prevalence_overall.csv, prevalence_by_group.csv, logistic_adjusted_or.csv, results.json
figures/            prevalence by BMI category, prevalence by age group
tests/              derivation logic, prevalence estimator, results.json consistency
REPORT.md           the written analysis
```

## Definitions

- **Hypertension**: mean SBP ≥ 130 mmHg **or** mean DBP ≥ 80 mmHg (2017 ACC/AHA)
  **or** currently taking antihypertensive medication (`BPQ050A`).
- **Obesity**: BMI ≥ 30 kg/m².
- **Sample**: adults aged ≥ 18 with a positive MEC weight and valid BP & BMI.

## Scope

A descriptive/reproducibility analysis on public data — it recovers
well-established associations (which is the point: the numbers match the
published literature, so the pipeline is trustworthy). It is **not** a novel
causal claim. Data: CDC/NCHS, public domain.

## License

MIT © Laura Piñero Roig
