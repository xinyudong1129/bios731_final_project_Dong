library(here)
library(dplyr)
library(readr)
library(lubridate)
library(data.table)
library(tibble)

source(here::here("R", "utils_paths.R"))
source(here::here("R", "preprocess_cgm.R"))
source(here::here("R", "feature_engineering.R"))

make_project_dirs()

raw_cgm <- read_cgm_data(here::here("data", "raw", "cgm.csv"))

analysis_data <- build_analysis_dataset(
  raw_data = raw_cgm,
  interval_mins = 5,
  window_size = 12,
  horizon_steps = 1
)

save_rds_here(analysis_data, "data", "processed", "analysis_data.rds")
write_csv_here(analysis_data, "data", "processed", "analysis_data.csv")