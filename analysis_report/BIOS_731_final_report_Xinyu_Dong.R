
source(here::here("scripts", "00_make_dirs.R"))
source(here::here("scripts", "01_build_analysis_data.R"))
source(here::here("scripts", "02_run_simulation_study.R"))
source(here::here("scripts", "03_run_real_data_analysis.R"))
source(here::here("scripts", "04_fit_bayes_real_full_data.R"))

readr::read_csv(here::here("results", "metrics", "real_data_cv_results.csv"))
readr::read_csv(here::here("results", "metrics", "real_data_cv_summary.csv"))

source(here::here("scripts", "05_summarize_results.R"))
