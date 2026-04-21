library(here)
library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(pROC)
library(PRROC)
library(rstanarm)

source(here::here("R", "utils_paths.R"))
source(here::here("R", "model_glm_em.R"))
source(here::here("R", "model_bayes.R"))
source(here::here("R", "evaluation.R"))
source(here::here("R", "pipeline_helpers.R"))

analysis_data <- read_rds_here("data", "processed", "analysis_data.rds")
folds <- make_group_folds(analysis_data, v = 5, seed = 123)

safe_eval <- function(expr, method_name, fold) {
  tryCatch(
    expr,
    error = function(e) {
      tibble::tibble(
        method = method_name,
        auroc = NA_real_,
        auprc = NA_real_,
        rmse_duration = NA_real_,
        coverage_95 = NA_real_,
        fold = fold,
        error_message = conditionMessage(e)
      )
    }
  )
}

cv_results <- purrr::map_dfr(seq_len(5), function(f) {
  sp <- split_by_fold(analysis_data, folds, fold = f)
  train <- sp$train
  test  <- sp$test
  
  if (length(unique(stats::na.omit(train$y_onset))) < 2) {
    return(tibble::tibble(
      method = c("glm_em", "bayes_mcmc"),
      auroc = NA_real_,
      auprc = NA_real_,
      rmse_duration = NA_real_,
      coverage_95 = NA_real_,
      fold = f,
      error_message = "Training fold has only one y_onset class."
    ))
  }
  
  res_glm <- safe_eval({
    fit_em_onset <- fit_glm_em_onset(train)
    fit_em_dur   <- fit_glm_em_duration(train)
    
    p_em <- predict_glm_em_onset(fit_em_onset, test)
    d_em <- predict_glm_em_duration(fit_em_dur, test)
    
    evaluate_predictions(
      data = test,
      onset_prob = p_em,
      duration_pred = d_em,
      method_name = "glm_em"
    ) %>%
      dplyr::mutate(fold = f, error_message = NA_character_)
  }, "glm_em", f)
  
  res_bayes <- safe_eval({
    # -------------------------
    # Bayesian onset model
    # -------------------------
    fit_mcmc_onset <- fit_bayes_mcmc_onset(
      train_data = train,
      chains = 2,
      iter = 1000,
      seed = 1000 + f
    )
    
    p_mcmc <- predict_bayes_onset(fit_mcmc_onset, test)
    
    # -------------------------
    # Bayesian duration model
    # -------------------------
    fit_mcmc_dur <- NULL
    d_mcmc <- rep(NA_real_, nrow(test))
    duration_lower <- NULL
    duration_upper <- NULL
    
    train_pos <- train %>%
      dplyr::filter(y_onset == 1, !is.na(y_duration))
    
    # only try duration model if enough positive rows
    if (nrow(train_pos) >= 30) {
      fit_mcmc_dur <- tryCatch(
        fit_bayes_mcmc_duration(
          train_data = train,
          chains = 2,
          iter = 1000,
          seed = 2000 + f
        ),
        error = function(e) NULL
      )
      
      if (!is.null(fit_mcmc_dur)) {
        d_mcmc <- tryCatch(
          predict_bayes_duration(fit_mcmc_dur, test),
          error = function(e) rep(NA_real_, nrow(test))
        )
        
        duration_pi <- tryCatch(
          posterior_interval_duration(fit_mcmc_dur, test),
          error = function(e) NULL
        )
        
        if (!is.null(duration_pi)) {
          duration_lower <- duration_pi$lower
          duration_upper <- duration_pi$upper
        }
      }
    }
    
    evaluate_predictions(
      data = test,
      onset_prob = p_mcmc,
      duration_pred = d_mcmc,
      duration_lower = duration_lower,
      duration_upper = duration_upper,
      method_name = "bayes_mcmc"
    ) %>%
      dplyr::mutate(fold = f, error_message = NA_character_)
  }, "bayes_mcmc", f)
  
  dplyr::bind_rows(res_glm, res_bayes)
})

write_csv_here(cv_results, "results", "metrics", "real_data_cv_results.csv")

cv_summary <- cv_results %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    auroc_mean = mean(auroc, na.rm = TRUE),
    auprc_mean = mean(auprc, na.rm = TRUE),
    rmse_mean = mean(rmse_duration, na.rm = TRUE),
    coverage_mean = mean(coverage_95, na.rm = TRUE),
    n_folds_with_auroc = sum(!is.na(auroc)),
    n_folds_with_auprc = sum(!is.na(auprc)),
    n_folds_with_rmse = sum(!is.na(rmse_duration)),
    .groups = "drop"
  )

write_csv_here(cv_summary, "results", "metrics", "real_data_cv_summary.csv")