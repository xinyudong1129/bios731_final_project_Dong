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
  window_size = 6,
  horizon_steps = 1,
  max_missing_in_window = 2,
  time_tolerance = 1
)

table(analysis_data$y_onset, useNA = "ifany")
mean(analysis_data$y_onset == 1, na.rm = TRUE)

save_rds_here(analysis_data, "data", "processed", "analysis_data.rds")
write_csv_here(analysis_data, "data", "processed", "analysis_data.csv")