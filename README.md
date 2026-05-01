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

###  Core Function Library (`R/`)

Contains reusable functions used throughout the project:

- `preprocess_cgm.R`  
  Cleans and formats raw CGM data (time alignment, handling gaps, standardization).

- `feature_engineering.R`  
  Constructs sliding-window features from CGM time series (mean, SD, slope, etc.).

- `model_glm_em.R`  
  Implements the baseline GLM (logistic regression) for onset prediction.

- `model_bayes.R`  
  Implements the Bayesian logistic regression model (MCMC-based).

- `evaluation.R`  
  Computes evaluation metrics (AUROC, AUPRC, RMSE).

- `pipeline_helpers.R`  
  Utility functions for orchestrating the pipeline (splits, wrappers, helpers).

- `simulate_cgm.R`  
  Generates synthetic CGM data for testing and simulation experiments.

- `utils_paths.R`  
  Handles file paths using `here::here()` for reproducibility across systems.

---

###  Pipeline Scripts (`scripts/`)

These scripts define the **end-to-end workflow**. They are designed to be run **in order**.

- `00_make_dirs.R`  
  Initializes the project structure by creating required folders (`results/`, `outputs/`, etc.).

- `01_build_analysis_data.R`  
  Loads raw CGM data, applies preprocessing, and generates the final analysis dataset with features and labels.

- `02_run_simulation_study.R` *(optional)*  
  Runs simulation experiments using synthetic CGM data to validate modeling approaches.

- `03_run_real_data_analysis.R`  
  Performs the main analysis on real CGM data:
  - Trains GLM and Bayesian models
  - Runs 5-fold grouped cross-validation
  - Saves predictions and metrics

- `04_fit_bayes_real_full_data.R`  
  Fits the Bayesian model on the **full dataset** (no CV) to produce final predicted probabilities used for risk stratification.

- `05_summarize_results.R`  
  Aggregates outputs into final tables and figures:
  - Cross-validation summaries
  - Risk decile tables
  - Plots used in the report

---

###  Data (`data/`)

- `raw/`  
  Original CGM data (unaltered)

- `processed/`  
  Cleaned and structured datasets ready for modeling

---

###  Results (`results/`)

Intermediate outputs generated during modeling:

- `metrics/`  
  Cross-validation performance results (AUROC, AUPRC, RMSE)

- `predictions/`  
  Model predictions (probabilities, durations)

- `models/`  
  Saved fitted model objects

- `simulation/`  
  Outputs from simulation experiments

---

###  Final Outputs (`outputs/`)

Materials used directly in the report:

- `figures/`  
  Plots (CV performance, risk deciles, duration scatter)

- `tables/`  
  Summary tables (CV results, deciles, top-risk windows)

---

###  Reports (`analysis_report/`)

- Final written report (R Markdown / PDF)
- Presentation slides

---

##  Workflow Summary

The full pipeline follows this order:

1. **Setup** → `00_make_dirs.R`  
2. **Data preparation** → `01_build_analysis_data.R`  
3. *(Optional)* Simulation → `02_run_simulation_study.R`  
4. **Model training + CV** → `03_run_real_data_analysis.R`  
5. **Final Bayesian model** → `04_fit_bayes_real_full_data.R`  
6. **Summarization + plots** → `05_summarize_results.R`

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




