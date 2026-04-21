library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)

source(here::here("R", "utils_paths.R"))

dir.create(here::here("outputs", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("outputs", "tables"), recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------
# 1. Read cross-validation results
# --------------------------------------------------
cv_results <- readr::read_csv(
  here::here("results", "metrics", "real_data_cv_results.csv"),
  show_col_types = FALSE
)

# More honest summary
cv_summary <- cv_results %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    auroc_mean = mean(auroc, na.rm = TRUE),
    auroc_sd   = sd(auroc, na.rm = TRUE),
    auprc_mean = mean(auprc, na.rm = TRUE),
    auprc_sd   = sd(auprc, na.rm = TRUE),
    rmse_mean  = mean(rmse_duration, na.rm = TRUE),
    rmse_sd    = sd(rmse_duration, na.rm = TRUE),
    n_folds_with_auroc = sum(!is.na(auroc)),
    n_folds_with_auprc = sum(!is.na(auprc)),
    n_folds_with_rmse  = sum(!is.na(rmse_duration)),
    .groups = "drop"
  )

readr::write_csv(
  cv_summary,
  here::here("outputs", "tables", "cv_summary_clean.csv")
)

# --------------------------------------------------
# 2. Pretty method labels
# --------------------------------------------------
cv_results <- cv_results %>%
  dplyr::mutate(
    method_label = dplyr::case_when(
      method == "glm_em" ~ "GLM baseline",
      method == "bayes_mcmc" ~ "Bayesian MCMC",
      method == "bayes_vb" ~ "Bayesian VB",
      TRUE ~ method
    )
  )

cv_summary <- cv_summary %>%
  dplyr::mutate(
    method_label = dplyr::case_when(
      method == "glm_em" ~ "GLM baseline",
      method == "bayes_mcmc" ~ "Bayesian MCMC",
      method == "bayes_vb" ~ "Bayesian VB",
      TRUE ~ method
    )
  )

# --------------------------------------------------
# 3. Bar plot: mean AUROC and AUPRC by method
# --------------------------------------------------
metric_bar_data <- cv_summary %>%
  dplyr::select(method_label, auroc_mean, auprc_mean) %>%
  tidyr::pivot_longer(
    cols = c(auroc_mean, auprc_mean),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           auroc_mean = "AUROC",
                           auprc_mean = "AUPRC"
    )
  )

p_bar <- ggplot(metric_bar_data, aes(x = method_label, y = value)) +
  geom_col() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Cross-validated onset prediction performance",
    x = NULL,
    y = "Mean metric value"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(
  filename = here::here("outputs", "figures", "cv_bar_metrics.png"),
  plot = p_bar,
  width = 8,
  height = 4.5,
  dpi = 300
)

# --------------------------------------------------
# 4. Fold-level plot for AUROC and AUPRC
# --------------------------------------------------
metric_fold_data <- cv_results %>%
  dplyr::select(method_label, fold, auroc, auprc) %>%
  tidyr::pivot_longer(
    cols = c(auroc, auprc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(metric,
                           auroc = "AUROC",
                           auprc = "AUPRC"
    )
  )

p_folds <- ggplot(metric_fold_data, aes(x = factor(fold), y = value, group = method_label)) +
  geom_point(size = 2) +
  geom_line(aes(color = method_label), linewidth = 0.7, na.rm = TRUE) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    title = "Fold-specific performance",
    x = "Fold",
    y = "Metric value",
    color = "Method"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = here::here("outputs", "figures", "cv_fold_metrics.png"),
  plot = p_folds,
  width = 8,
  height = 4.5,
  dpi = 300
)

# --------------------------------------------------
# 5. RMSE plot if available
# --------------------------------------------------
rmse_data <- cv_results %>%
  dplyr::filter(!is.na(rmse_duration)) %>%
  dplyr::mutate(
    method_label = dplyr::case_when(
      method == "glm_em" ~ "GLM baseline",
      method == "bayes_mcmc" ~ "Bayesian MCMC",
      method == "bayes_vb" ~ "Bayesian VB",
      TRUE ~ method
    )
  )

if (nrow(rmse_data) > 0) {
  p_rmse <- ggplot(rmse_data, aes(x = factor(fold), y = rmse_duration, color = method_label, group = method_label)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    labs(
      title = "Fold-specific duration RMSE",
      x = "Fold",
      y = "RMSE",
      color = "Method"
    ) +
    theme_minimal(base_size = 13)
  
  ggsave(
    filename = here::here("outputs", "figures", "cv_rmse.png"),
    plot = p_rmse,
    width = 7,
    height = 4.5,
    dpi = 300
  )
}

# --------------------------------------------------
# 6. Full-data Bayesian onset predictions
# --------------------------------------------------
onset_pred_path <- here::here("results", "predictions", "bayes_mcmc_onset_full_data_predictions.csv")

if (file.exists(onset_pred_path)) {
  onset_pred <- readr::read_csv(onset_pred_path, show_col_types = FALSE)
  
  # Histogram of predicted onset probabilities
  p_hist <- ggplot(onset_pred, aes(x = onset_prob)) +
    geom_histogram(bins = 40) +
    labs(
      title = "Distribution of Bayesian predicted onset probabilities",
      x = "Predicted probability of missing interval onset",
      y = "Count"
    ) +
    theme_minimal(base_size = 13)
  
  ggsave(
    filename = here::here("outputs", "figures", "bayes_onset_prob_hist.png"),
    plot = p_hist,
    width = 7,
    height = 4.5,
    dpi = 300
  )
  
  # Top-risk windows
  top_risk <- onset_pred %>%
    dplyr::arrange(dplyr::desc(onset_prob)) %>%
    dplyr::slice_head(n = 20)
  
  readr::write_csv(
    top_risk,
    here::here("outputs", "tables", "top20_bayesian_onset_risk.csv")
  )
  
  # Decile enrichment
  onset_deciles <- onset_pred %>%
    dplyr::mutate(
      risk_decile = dplyr::ntile(onset_prob, 10)
    ) %>%
    dplyr::group_by(risk_decile) %>%
    dplyr::summarise(
      n = dplyr::n(),
      observed_onset_rate = mean(y_onset == 1, na.rm = TRUE),
      mean_predicted_prob = mean(onset_prob, na.rm = TRUE),
      .groups = "drop"
    )
  
  readr::write_csv(
    onset_deciles,
    here::here("outputs", "tables", "bayes_onset_deciles.csv")
  )
  
  p_decile <- ggplot(onset_deciles, aes(x = factor(risk_decile), y = observed_onset_rate)) +
    geom_col() +
    labs(
      title = "Observed onset rate by predicted-risk decile",
      x = "Predicted risk decile",
      y = "Observed onset rate"
    ) +
    theme_minimal(base_size = 13)
  
  ggsave(
    filename = here::here("outputs", "figures", "bayes_onset_decile_plot.png"),
    plot = p_decile,
    width = 7,
    height = 4.5,
    dpi = 300
  )
}

# --------------------------------------------------
# 7. Full-data Bayesian duration predictions
# --------------------------------------------------
duration_pred_path <- here::here("results", "predictions", "bayes_mcmc_duration_full_data_predictions.csv")

if (file.exists(duration_pred_path)) {
  duration_pred <- readr::read_csv(duration_pred_path, show_col_types = FALSE)
  
  duration_eval <- duration_pred %>%
    dplyr::filter(y_onset == 1, !is.na(y_duration), !is.na(duration_pred))
  
  if (nrow(duration_eval) > 0) {
    p_scatter <- ggplot(duration_eval, aes(x = y_duration, y = duration_pred)) +
      geom_point(alpha = 0.7) +
      geom_abline(intercept = 0, slope = 1, linetype = 2) +
      labs(
        title = "Observed vs predicted missing duration",
        x = "Observed duration",
        y = "Predicted duration"
      ) +
      theme_minimal(base_size = 13)
    
    ggsave(
      filename = here::here("outputs", "figures", "bayes_duration_scatter.png"),
      plot = p_scatter,
      width = 6.5,
      height = 5,
      dpi = 300
    )
  }
}