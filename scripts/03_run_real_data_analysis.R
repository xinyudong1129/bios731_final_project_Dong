library(here)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(pROC)
library(PRROC)
library(flexmix)
library(rstanarm)

source(here::here("R", "utils_paths.R"))
source(here::here("R", "model_glm_em.R"))
source(here::here("R", "model_bayes.R"))
source(here::here("R", "evaluation.R"))
source(here::here("R", "pipeline_helpers.R"))

analysis_data <- read_rds_here("data", "processed", "analysis_data.rds")
folds <- make_group_folds(analysis_data, v = 5, seed = 123)

cv_results <- purrr::map_dfr(seq_len(5), function(f) {
  sp <- split_by_fold(analysis_data, folds, fold = f)
  train <- sp$train
  test  <- sp$test
  
  fit_em_onset <- fit_glm_em_onset(train)
  fit_em_dur   <- fit_glm_em_duration(train)
  p_em <- predict_glm_em_onset(fit_em_onset, test)
  d_em <- predict_glm_em_duration(fit_em_dur, test)
  res_em <- evaluate_predictions(test, p_em, d_em, method_name = "glm_em")
  
  fit_mcmc_onset <- fit_bayes_mcmc_onset(train, chains = 2, iter = 1000, seed = 11 + f)
  fit_mcmc_dur   <- fit_bayes_mcmc_duration(train, chains = 2, iter = 1000, seed = 21 + f)
  p_mcmc <- predict_bayes_onset(fit_mcmc_onset, test)
  d_mcmc <- predict_bayes_duration(fit_mcmc_dur, test)
  pi_mcmc <- posterior_interval_duration(fit_mcmc_dur, test)
  res_mcmc <- evaluate_predictions(
    test, p_mcmc, d_mcmc,
    duration_lower = pi_mcmc$lower,
    duration_upper = pi_mcmc$upper,
    method_name = "bayes_mcmc"
  )
  
  fit_vb_onset <- fit_bayes_vb_onset(train, seed = 31 + f)
  fit_vb_dur   <- fit_bayes_vb_duration(train, seed = 41 + f)
  p_vb <- predict_bayes_onset(fit_vb_onset, test)
  d_vb <- predict_bayes_duration(fit_vb_dur, test)
  pi_vb <- posterior_interval_duration(fit_vb_dur, test)
  res_vb <- evaluate_predictions(
    test, p_vb, d_vb,
    duration_lower = pi_vb$lower,
    duration_upper = pi_vb$upper,
    method_name = "bayes_vb"
  )
  
  dplyr::bind_rows(res_em, res_mcmc, res_vb) |>
    dplyr::mutate(fold = f)
})

write_csv_here(cv_results, "results", "metrics", "real_data_cv_results.csv")

cv_summary <- cv_results |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    auroc_mean = mean(auroc, na.rm = TRUE),
    auprc_mean = mean(auprc, na.rm = TRUE),
    rmse_mean = mean(rmse_duration, na.rm = TRUE),
    coverage_mean = mean(coverage_95, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_here(cv_summary, "results", "metrics", "real_data_cv_summary.csv")