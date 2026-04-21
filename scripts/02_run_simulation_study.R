library(here)
library(dplyr)
library(purrr)
library(tibble)
library(readr)

library(pROC)
library(PRROC)
library(flexmix)
library(rstanarm)
library(lubridate)
library(data.table)

source(here::here("R", "utils_paths.R"))
source(here::here("R", "preprocess_cgm.R"))
source(here::here("R", "feature_engineering.R"))
source(here::here("R", "simulate_cgm.R"))
source(here::here("R", "model_glm_em.R"))
source(here::here("R", "model_bayes.R"))
source(here::here("R", "evaluation.R"))
source(here::here("R", "pipeline_helpers.R"))

#plan(multisession)

run_one_simulation <- function(rep_id, n_patients = 50) {
  sim_raw <- simulate_cgm_panel(n_patients = n_patients, seed = 1000 + rep_id)
  
  sim_dat <- build_analysis_dataset(
    raw_data = sim_raw,
    interval_mins = 5,
    window_size = 12,
    horizon_steps = 1
  )
  
  folds <- make_group_folds(sim_dat, v = 5, seed = 2000 + rep_id)
  
  fold_results <- purrr::map_dfr(seq_len(5), function(f) {
    sp <- split_by_fold(sim_dat, folds, fold = f)
    train <- sp$train
    test  <- sp$test
    
    # skip unusable folds
    if (length(unique(train$y_onset)) < 2) {
      return(tibble::tibble(
        method = c("glm_em", "bayes_mcmc"),
        auroc = NA_real_,
        auprc = NA_real_,
        rmse_duration = NA_real_,
        coverage_95 = NA_real_,
        rep_id = rep_id,
        fold = f
      ))
    }
    
    safe_eval <- function(expr, method_name) {
      tryCatch(
        expr,
        error = function(e) {
          tibble::tibble(
            method = method_name,
            auroc = NA_real_,
            auprc = NA_real_,
            rmse_duration = NA_real_,
            coverage_95 = NA_real_,
            rep_id = rep_id,
            fold = f
          )
        }
      )
    }
    
    res_em <- safe_eval({
      fit_em_onset <- fit_glm_em_onset(train)
      fit_em_dur   <- fit_glm_em_duration(train)
      test <- test[, names(train), drop = FALSE]
      p_em <- predict_glm_em_onset(fit_em_onset, test)
      d_em <- predict_glm_em_duration(fit_em_dur, test)
      evaluate_predictions(test, p_em, d_em, method_name = "glm_em") %>%
        dplyr::mutate(rep_id = rep_id, fold = f)
    }, "glm_em")
    
    res_mcmc <- safe_eval({
      fit_mcmc_onset <- fit_bayes_mcmc_onset(train, chains = 4, iter = 3000, seed = 3000 + f)
      fit_mcmc_dur   <- fit_bayes_mcmc_duration(train, chains = 4, iter = 3000, seed = 4000 + f)
      test <- test[, names(train), drop = FALSE]
      p_mcmc <- predict_bayes_onset(fit_mcmc_onset, test)
      d_mcmc <- predict_bayes_duration(fit_mcmc_dur, test)
      pi_mcmc <- posterior_interval_duration(fit_mcmc_dur, test)
      evaluate_predictions(
        test, p_mcmc, d_mcmc,
        duration_lower = pi_mcmc$lower,
        duration_upper = pi_mcmc$upper,
        method_name = "bayes_mcmc"
      ) %>% dplyr::mutate(rep_id = rep_id, fold = f)
    }, "bayes_mcmc")
   
    dplyr::bind_rows(res_em, res_mcmc)
  })
  
  fold_results
}

sim_results <- purrr::map_dfr(
  .x = 1:5,
  .f = run_one_simulation
)

write_csv_here(sim_results, "results", "simulation", "simulation_results.csv")

sim_summary <- sim_results |>
  dplyr::group_by(method) |>
  dplyr::summarise(
    auroc_mean = mean(auroc, na.rm = TRUE),
    auprc_mean = mean(auprc, na.rm = TRUE),
    rmse_mean = mean(rmse_duration, na.rm = TRUE),
    coverage_mean = mean(coverage_95, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_here(sim_summary, "results", "simulation", "simulation_summary.csv")