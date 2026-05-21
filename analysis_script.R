# -------------------------------------------------------
# analyze_results.R
# Loads all simulation results and computes bias and variance
# -------------------------------------------------------

library(dplyr)
library(ggplot2)

# -------------------------------------------------------
# Parameter mapping
# -------------------------------------------------------
param_map <- data.frame(
  param_num  = 17:24,
  magnitude  = c(0, 0, 0, 0, 1, 1, 1, 1),
  noise      = c(0, 0, 1, 1, 0, 0, 1, 1),
  confounding = c(0, 1, 0, 1, 0, 1, 0, 1)
)

# -------------------------------------------------------
# Load all results
# -------------------------------------------------------
results_dir <- "./results"
files <- list.files(results_dir, pattern = "^results_run02.*\\.RData$", full.names = TRUE)

cat(sprintf("Loading %d result files...\n", length(files)))

all_results <- lapply(files, function(f) {
  e <- new.env()
  load(f, envir = e)
  e$results_df
})

results_df <- do.call(rbind, all_results)

# Merge with param map
results_df <- results_df %>%
  left_join(param_map, by = "param_num")

cat(sprintf("Total rows loaded: %d\n", nrow(results_df)))

# -------------------------------------------------------
# Summary: bias and variance by estimator and param
# -------------------------------------------------------
summary_df <- results_df %>%
  group_by(param_num, magnitude, noise, confounding, estimator, sigma_type) %>%
  summarise(
    mean_bias     = mean(bias, na.rm = TRUE),
    abs_bias      = mean(abs(bias), na.rm = TRUE),
    variance      = var(ATE_estimate, na.rm = TRUE),
    rmse          = sqrt(mean(bias^2, na.rm = TRUE)),
    n_sims        = n(),
    .groups = "drop"
  )

print(summary_df)

# -------------------------------------------------------
# Bias reduction: AIPW vs HOIF by confounding
# -------------------------------------------------------
bias_by_confounding <- results_df %>%
  group_by(confounding, estimator) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    abs_bias  = mean(abs(bias), na.rm = TRUE),
    .groups   = "drop"
  )

# -------------------------------------------------------
# Bias reduction: AIPW vs HOIF by noise
# -------------------------------------------------------
bias_by_noise <- results_df %>%
  group_by(noise, estimator) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    abs_bias  = mean(abs(bias), na.rm = TRUE),
    .groups   = "drop"
  )

# -------------------------------------------------------
# Bias reduction: AIPW vs HOIF by magnitude
# -------------------------------------------------------
bias_by_magnitude <- results_df %>%
  group_by(magnitude, estimator) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    abs_bias  = mean(abs(bias), na.rm = TRUE),
    .groups   = "drop"
  )

# -------------------------------------------------------
# Save all
# -------------------------------------------------------
save(summary_df, results_df, bias_by_confounding, bias_by_noise, bias_by_magnitude,
     file = "./results/analysis_results_run02.RData")
cat("\nSaved to analysis_results_run02.RData\n")

