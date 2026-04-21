library(here)
library(dplyr)
library(readr)
library(rstanarm)
library(tibble)

source(here::here("R", "utils_paths.R"))
source(here::here("R", "model_bayes.R"))

analysis_data <- read_rds_here("data", "processed", "analysis_data.rds")

# --------------------------------------------------
# Keep only columns needed by the Bayesian models
# --------------------------------------------------
onset_vars <- c(
  "patient_id", "end_time", "y_onset",
  "mean_glucose", "sd_glucose", "slope_glucose",
  "delta_last_first", "pct_high", "pct_low"
)

duration_vars <- c(
  "patient_id", "end_time", "y_onset", "y_duration",
  "mean_glucose", "sd_glucose", "slope_glucose",
  "delta_last_first", "pct_high", "pct_low"
)

# --------------------------------------------------
# Clean onset dataset
# --------------------------------------------------
analysis_onset <- analysis_data %>%
  dplyr::select(dplyr::all_of(onset_vars)) %>%
  dplyr::filter(
    !is.na(y_onset),
    is.finite(mean_glucose),
    is.finite(sd_glucose),
    is.finite(slope_glucose),
    is.finite(delta_last_first),
    is.finite(pct_high),
    is.finite(pct_low)
  )

cat("Rows in raw analysis_data:", nrow(analysis_data), "\n")
cat("Rows in cleaned onset data:", nrow(analysis_onset), "\n")
cat("Onset class counts:\n")
print(table(analysis_onset$y_onset, useNA = "ifany"))

if (nrow(analysis_onset) == 0) {
  stop("No usable rows remain for Bayesian onset model after cleaning.")
}

if (length(unique(analysis_onset$y_onset)) < 2) {
  stop("Bayesian onset model requires both y_onset classes after cleaning.")
}

# --------------------------------------------------
# Fit Bayesian onset model
# --------------------------------------------------
cat("\nFitting Bayesian onset model...\n")

fit_onset <- fit_bayes_mcmc_onset(
  train_data = analysis_onset,
  chains = 4,
  iter = 2000,
  seed = 123
)

save_rds_here(
  fit_onset,
  "results", "models", "bayes_mcmc_onset_full_data.rds"
)

onset_prob <- predict_bayes_onset(fit_onset, analysis_onset)

onset_predictions <- analysis_onset %>%
  dplyr::select(patient_id, end_time, y_onset) %>%
  dplyr::mutate(onset_prob = onset_prob)

write_csv_here(
  onset_predictions,
  "results", "predictions", "bayes_mcmc_onset_full_data_predictions.csv"
)

# --------------------------------------------------
# Clean duration dataset
# --------------------------------------------------
analysis_duration <- analysis_data %>%
  dplyr::select(dplyr::all_of(duration_vars)) %>%
  dplyr::filter(
    y_onset == 1,
    !is.na(y_duration),
    is.finite(y_duration),
    is.finite(mean_glucose),
    is.finite(sd_glucose),
    is.finite(slope_glucose),
    is.finite(delta_last_first),
    is.finite(pct_high),
    is.finite(pct_low)
  )

cat("\nRows in cleaned duration data:", nrow(analysis_duration), "\n")

if (nrow(analysis_duration) < 10) {
  message("Not enough positive-duration rows to fit full-data Bayesian duration model.")
} else {
  cat("Fitting Bayesian duration model...\n")
  
  fit_duration <- fit_bayes_mcmc_duration(
    train_data = analysis_duration,
    chains = 4,
    iter = 2000,
    seed = 456
  )
  
  if (is.null(fit_duration)) {
    message("Duration fit returned NULL.")
  } else {
    save_rds_here(
      fit_duration,
      "results", "models", "bayes_mcmc_duration_full_data.rds"
    )
    
    duration_pred <- predict_bayes_duration(fit_duration, analysis_duration)
    duration_pi <- posterior_interval_duration(fit_duration, analysis_duration)
    
    duration_predictions <- analysis_duration %>%
      dplyr::select(patient_id, end_time, y_onset, y_duration) %>%
      dplyr::mutate(
        duration_pred = duration_pred,
        duration_lower = duration_pi$lower,
        duration_upper = duration_pi$upper
      )
    
    write_csv_here(
      duration_predictions,
      "results", "predictions", "bayes_mcmc_duration_full_data_predictions.csv"
    )
  }
}