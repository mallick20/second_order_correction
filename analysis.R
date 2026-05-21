library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)


# df <- read.csv("results_for_r.csv")
load("analysis_results_run02.RData")
df <- results_df
df <- df %>% mutate(abs_bias = abs(bias))

### Failure check

failure_check <- df %>%
  mutate(
    abs_bias = abs(ATE_estimate - true_ATE),
    numerical_failure = abs_bias > 1e4 | is.infinite(abs_bias) | is.na(abs_bias)
  ) %>%
  group_by(sigma_type) %>%
  summarise(
    n = n(),
    failures = sum(numerical_failure),
    failure_rate = mean(numerical_failure),
    max_abs_bias = max(abs_bias, na.rm = TRUE),
    p99_abs_bias = quantile(abs_bias, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

print(failure_check)



### Plot checking distribution of AIPW and abs_biases across runs

plot_df_1 <- df %>%
  filter(
    is.na(sigma_type) |
      str_detect(sigma_type, "^tr")
  ) %>%
  mutate(
    estimator = case_when(
      is.na(sigma_type) ~ "AIPW",
      sigma_type == "tr1" ~ "HOIF d1",
      sigma_type == "tr2" ~ "HOIF d2",
      sigma_type == "tr3" ~ "HOIF d3"
    )
  )

ggplot(
  plot_df_1,
  aes(
    x = estimator,
    y = abs_bias,
    fill = estimator
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  labs(
    title = "Log Absolute Bias Across IID Simulation Runs",
    x = "Estimator",
    y = "|Estimated ATE - True ATE| (log scale)"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))

plot_df_1 %>% filter(estimator %in%  c("AIPW","HOIF d1")) %>%  ggplot(
  aes(
    x = estimator,
    y = bias,
    fill = estimator
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    title = "Bias for AIPW and HOIF d1",
    x = "Estimator",
    y = "Estimated ATE - True ATE"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))


library(dplyr)

summary_table <- results_df %>%
  filter(
    is.na(sigma_type) |
      str_detect(sigma_type, "^tr")
  ) %>%
  group_by(estimator) %>%
  summarise(
    bias = mean(ATE_estimate - true_ATE),
    abs_bias = mean(abs(ATE_estimate - true_ATE)),
    variance = var(ATE_estimate - true_ATE),
    sd = sd(ATE_estimate - true_ATE),
    rmse = sqrt(mean((ATE_estimate - true_ATE)^2))
  )

print(summary_table)


#### Plot distribution of estimations across confounding and non-confounding
confounding_plot_df <- df %>%
  filter(
    is.na(sigma_type) |
      str_detect(sigma_type, "^tr")
  ) %>%
  mutate(
    estimator = case_when(
      is.na(sigma_type) ~ "AIPW",
      sigma_type == "tr1" ~ "HOIF d1",
      sigma_type == "tr2" ~ "HOIF d2",
      sigma_type == "tr3" ~ "HOIF d3"
    ),
    confounding_label = case_when(
      confounding == 0 ~ "Weak confounding",
      confounding == 1 ~ "Strong confounding"
    )
  )

ggplot(
  confounding_plot_df,
  aes(
    x = estimator,
    y = abs_bias,
    fill = estimator
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  facet_wrap(~ confounding_label) +
  labs(
    title = "Absolute Bias by Confounding Strength",
    x = "Estimator",
    y = "Absolute Bias (log scale)"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))



#### Plot distribution of estimations across lower magnitude effect and higher magnitude effect
magnitude_plot_df <- df %>%
  filter(
    is.na(sigma_type) |
      str_detect(sigma_type, "^tr")
  ) %>%
  mutate(
    estimator = case_when(
      is.na(sigma_type) ~ "AIPW",
      sigma_type == "tr1" ~ "HOIF d1",
      sigma_type == "tr2" ~ "HOIF d2",
      sigma_type == "tr3" ~ "HOIF d3"
    ),
    magnitude_label = case_when(
      magnitude == 0 ~ "Low effect magnitude",
      magnitude == 1 ~ "High effect magnitude"
    )
  )

ggplot(
  magnitude_plot_df,
  aes(
    x = estimator,
    y = abs_bias,
    fill = estimator
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  facet_wrap(~ magnitude_label) +
  labs(
    title = "Absolute Bias by Effect Magnitude",
    x = "Estimator",
    y = "Absolute Bias (log scale)"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))


#### Plot 
noise_plot_df <- df %>%
  filter(
    is.na(sigma_type) |
      str_detect(sigma_type, "^tr")
  ) %>%
  mutate(
    estimator = case_when(
      is.na(sigma_type) ~ "AIPW",
      sigma_type == "tr1" ~ "HOIF d1",
      sigma_type == "tr2" ~ "HOIF d2",
      sigma_type == "tr3" ~ "HOIF d3"
    ),
    noise_label = case_when(
      noise == 0 ~ "Low noise",
      noise == 1 ~ "High noise"
    )
  )

ggplot(
  noise_plot_df,
  aes(
    x = estimator,
    y = abs_bias,
    fill = estimator
  )
) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_y_log10() +
  facet_wrap(~ noise_label) +
  labs(
    title = "Absolute Bias by Noise Level",
    x = "Estimator",
    y = "Absolute Bias (log scale)"
  ) +
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5))



#### Overall difference
config_plot_df <- df %>%
  filter(
    is.na(sigma_type) |
      str_detect(sigma_type, "^tr")
  ) %>%
  mutate(
    estimator = case_when(
      is.na(sigma_type) ~ "AIPW",
      sigma_type == "tr1" ~ "HOIF d1",
      sigma_type == "tr2" ~ "HOIF d2",
      sigma_type == "tr3" ~ "HOIF d3"
    ),
    config_label = paste0(
      "P", param_num,
      "\nM=", magnitude,
      ", N=", noise,
      ", C=", confounding
    )
  )


ggplot(
  config_plot_df,
  aes(
    x = config_label,
    y = abs_bias,
    fill = estimator
  )
) +
  geom_boxplot(outlier.alpha = 0.25) +
  scale_y_log10() +
  labs(
    title = "Absolute Error by IID Configuration",
    x = "Configuration",
    y = "Absolute Error (log scale)",
    fill = "Estimator"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# Keep only AIPW and standard covariance HOIF
win_df <- df %>%
  filter(is.na(sigma_type) | str_detect(sigma_type, "^tr")) %>%
  mutate(
    estimator = case_when(
      is.na(sigma_type) ~ "AIPW",
      sigma_type == "tr1" ~ "HOIF d1",
      sigma_type == "tr2" ~ "HOIF d2",
      sigma_type == "tr3" ~ "HOIF d3"
    )
  ) %>%
  select(param_num, sim_num, magnitude, noise, confounding, estimator, abs_bias)

# Wide format: one row per simulation/config
win_wide <- win_df %>%
  pivot_wider(
    names_from = estimator,
    values_from = abs_bias
  )

# Win rate: HOIF abs_error < AIPW abs_error
win_summary <- win_wide %>%
  summarise(
    across(
      starts_with("HOIF"),
      ~ mean(.x < AIPW, na.rm = TRUE),
      .names = "{.col}"
    ),
    .by = c(param_num, magnitude, noise, confounding)
  ) %>%
  pivot_longer(
    cols = starts_with("HOIF"),
    names_to = "estimator",
    values_to = "win_rate"
  ) %>%
  mutate(
    config_label = paste0(
      "P", param_num,
      "\nM=", magnitude,
      ", N=", noise,
      ", C=", confounding
    )
  )

ggplot(
  win_summary,
  aes(
    x = config_label,
    y = win_rate,
    fill = estimator
  )
) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(
    title = "HOIF Win Rate Relative to AIPW",
    x = "Configuration",
    y = "Proportion of runs where HOIF absolute error < AIPW",
    fill = "Estimator"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


