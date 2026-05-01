#  bios731_final_project_Dong

##  Predicting Missing Intervals in Continuous Glucose Monitoring (CGM) Data

This repository contains a full, reproducible analysis pipeline for predicting **missing interval onset and duration** in CGM data.
---

##  Overview

Continuous Glucose Monitoring (CGM) data are high-frequency time series widely used in diabetes research. However, **missing intervals** frequently occur due to sensor failure, detachment, or calibration issues.

This project addresses two key prediction tasks:

1. **Onset Prediction (Binary)**
   - Predict whether a missing interval begins after a given CGM window

2. **Duration Prediction (Continuous)**
   - Predict how long the missing interval lasts

---

##  Methodology

###  Feature Engineering
- Sliding windows over CGM time series
- Extract features capturing:
  - Glucose level (mean)
  - Variability (SD, CV)
  - Trend (slope)
  - Extremes (min/max)
  - Short-term dynamics

###  Models

**Onset model (rare-event classification):**
- GLM baseline (logistic regression)
- Bayesian logistic regression (MCMC)

**Duration model (conditional regression):**
- Linear regression on onset-positive windows

###  Evaluation
- 5-fold grouped cross-validation (by patient)
- Metrics:
  - AUROC
  - AUPRC (critical for rare events)
  - RMSE (duration)

---

##  Repository Structure

├── R/ # Core reusable R functions
│ ├── feature_engineering.R
│ ├── model_bayes.R
│ ├── model_glm_em.R
│ ├── preprocess_cgm.R
│ ├── evaluation.R
│ ├── pipeline_helpers.R
│ ├── simulate_cgm.R
│ └── utils_paths.R
│
├── scripts/ # Pipeline execution scripts (ordered)
│ ├── 00_make_dirs.R
│ ├── 01_build_analysis_data.R
│ ├── 02_run_simulation_study.R
│ ├── 03_run_real_data_analysis.R
│ ├── 04_fit_bayes_real_full_data.R
│ └── 05_summarize_results.R
│
├── data/
│ ├── raw/ # Original CGM data
│ └── processed/ # Cleaned data
│
├── results/
│ ├── metrics/ # Cross-validation results
│ ├── predictions/ # Model predictions
│ ├── models/ # Saved models
│ └── simulation/ # Simulation outputs
│
├── outputs/
│ ├── figures/ # Plots used in report
│ └── tables/ # Summary tables
│
├── analysis_report/ # Final report + slides
├── BIOS 731 final project.Rproj
└── .gitignore


## For reproduction:

Run scripts in order:

```r
source("scripts/00_make_dirs.R")
source("scripts/01_build_analysis_data.R")
source("scripts/02_run_simulation_study.R")   # optional
source("scripts/03_run_real_data_analysis.R")
source("scripts/04_fit_bayes_real_full_data.R")
source("scripts/05_summarize_results.R")
```

## Outputs
`outputs/figures/`: plots (CV performance, deciles, duration scatter)
`outputs/tables/`: summaries and results tables




