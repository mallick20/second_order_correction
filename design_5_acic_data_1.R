source('./src/dependencies_min.R')
source('./src/estimation_functions.R')


# Loading data from acic
load("data/input_2017.RData")
load("data/parameters_2017.RData")

# Run once in your R console to install:
#   devtools::install_github("vdorie/aciccomp/2017", upgrade = "never")
library(aciccomp2017)

# -------------------------------------------------------
# Design Parameters
# -------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
param_num <- if (length(args) >= 1) as.integer(args[1]) else 1
sim_num   <- if (length(args) >= 2) as.integer(args[2]) else 1
run_id <- as.integer(args[3])
cat(sprintf("Running: param_num = %d | sim_num = %d\n", param_num, sim_num))
n_folds   <- 5   # number of cross-fitting folds
holdout_frac <- 0

degree_of_interactions <- c(1,2,3)
# Covariates used in the DGP (paper Section 3)
# Continuous: x1, x3, x43
# Binary:     x10, x14, x15, x24
# Categorical: x21
dgp_covariates  <- c("x_1", "x_3", "x_10", "x_14", "x_15", "x_21", "x_24", "x_43")
categorical_vars     <- c("x_10", "x_14", "x_15", "x_24")
continuous_vars <- c("x_1", "x_3", "x_43", "x_21")   # x21 treated as continuous (16 levels)

set.seed(123)

# -------------------------------------------------------
# Load data from aciccomp2017
#   dgp_2017() returns a plain data.frame with:
#   - covariates from input_2017 (n = 4,302 rows)
#   - z: treatment indicator
#   - y: outcome
# -------------------------------------------------------
data <- dgp_2017(param_num, sim_num) %>%
  rename(A = z, Y = y) 

full_data <- input_2017 %>% dplyr::select(all_of(dgp_covariates))
full_data <- cbind(full_data, data)

# Changing factors to numeric columns
factor_cols <- names(full_data)[sapply(full_data, is.factor)]
full_data[factor_cols] <- lapply(full_data[factor_cols], as.numeric)



# -------------------------------------------------------
# True ATE: tau(x) = xi * (x3*x24 + (x14-1) - (x15-1))
# xi lookup from paper Table 2
# -------------------------------------------------------
true_ATE <- mean(full_data$alpha)

cat(sprintf("True ATE (full data): %.4f\n", true_ATE))

# -------------------------------------------------------
# Split: 10% holdout, 90% for cross-fitting
# -------------------------------------------------------
n_total   <- nrow(full_data)
n_holdout <- floor(holdout_frac * n_total)

holdout_idx   <- sample(n_total, n_holdout, replace = FALSE)
crossfit_idx  <- setdiff(seq_len(n_total), holdout_idx)

data_holdout  <- full_data[holdout_idx, ]
data_crossfit <- full_data[crossfit_idx, ]
n_crossfit    <- nrow(data_crossfit)

cat(sprintf("Total rows: %d | Cross-fit: %d | Holdout: %d\n",
            n_total, nrow(data_crossfit), nrow(data_holdout)))

# -------------------------------------------------------
# 5-fold cross-fitting IDs on the 90% data
# -------------------------------------------------------
fold_ids <- get_CV_ids(data_crossfit, k = n_folds)

cat(sprintf("Fold sizes: %s\n",
            paste(table(fold_ids), collapse = " / ")))



# -------------------------------------------------------
# Nuisance model registry
# Add new models here as named entries
# -------------------------------------------------------
nuisance_models <- list(
  glm = list(
    params_trt = find_params_glm(
      binary_vars                = categorical_vars,
      continuous_vars            = continuous_vars,
      continuous_var_spline_knots = 2),
    params_outcome = find_params_glm(
      binary_vars                = c(categorical_vars, "A"),
      continuous_vars            = continuous_vars,
      continuous_var_spline_knots = 2),
    model_name = "glm"
  )
  # Add future models here, e.g.:
  # lasso = list(params_trt = ..., params_outcome = ..., model_name = "lasso")
  # random_forest = list(...)
)

# -------------------------------------------------------
# Cross-fitting loop over all models
# -------------------------------------------------------
cf_predictions <- lapply(nuisance_models, function(model) {
  prob_trt1          <- numeric(n_crossfit)
  prob_outcome1_trt1 <- numeric(n_crossfit)
  prob_outcome1_trt0 <- numeric(n_crossfit)
  
  for (k in 1:n_folds) {
    fold_est_idx   <- which(fold_ids == k)
    fold_train_idx <- which(fold_ids != k)
    fold_train     <- data_crossfit[fold_train_idx, ]
    fold_est       <- data_crossfit[fold_est_idx, ]
    
    prob_trt1[fold_est_idx] <- estimate_prob_individual_model(
      covariates_df = fold_train %>% dplyr::select(-c(Y, A)),
      label_vector  = fold_train$A,
      params        = model$params_trt,
      model         = model$model_name,
      predict_data  = fold_est %>% dplyr::select(-c(Y, A)))
    
    prob_outcome1_trt1[fold_est_idx] <- estimate_prob_individual_model(
      covariates_df = fold_train %>% dplyr::select(-Y),
      label_vector  = fold_train$Y,
      params        = model$params_outcome,
      model         = model$model_name,
      predict_data  = fold_est %>% dplyr::select(-Y) %>% mutate(A = 1))
    
    prob_outcome1_trt0[fold_est_idx] <- estimate_prob_individual_model(
      covariates_df = fold_train %>% dplyr::select(-Y),
      label_vector  = fold_train$Y,
      params        = model$params_outcome,
      model         = model$model_name,
      predict_data  = fold_est %>% dplyr::select(-Y) %>% mutate(A = 0))
  }
  
  list(prob_trt1          = prob_trt1,
       prob_outcome1_trt1 = prob_outcome1_trt1,
       prob_outcome1_trt0 = prob_outcome1_trt0)
})



# AIPW on cross-fitted predictions
aipw_results <- lapply(names(cf_predictions), function(model_name) {
  preds <- cf_predictions[[model_name]]
  
  list(
    estimate = estimate_AIPW(
      trt                = data_crossfit$A,
      outcome            = data_crossfit$Y,
      prob_trt1          = preds$prob_trt1,
      prob_outcome1_trt1 = preds$prob_outcome1_trt1,
      prob_outcome1_trt0 = preds$prob_outcome1_trt0),
    variance = estimate_var_AIPW(
      trt                = data_crossfit$A,
      outcome            = data_crossfit$Y,
      prob_trt1          = preds$prob_trt1,
      prob_outcome1_trt1 = preds$prob_outcome1_trt1,
      prob_outcome1_trt0 = preds$prob_outcome1_trt0)
  )
})
names(aipw_results) <- names(cf_predictions)




#### Second Order Estimator


#### Basis Creation - Type 1

basis_cf <- lapply(degree_of_interactions, function(k) {
  create_b_spline_basis(
    data                   = data_crossfit,
    binary_vars            = categorical_vars,
    continuous_vars        = continuous_vars,
    degree_of_interactions = k
  )
})


#### Sigma Creation
sigma_cf_eff1 <- c(
  'tr' = compute_sigma(basis = basis_cf, trt = data_crossfit$A),
  'nlshrink' = lapply(1:length(degree_of_interactions), function(i) {
    nlshrink_cov(data_crossfit$A * basis_cf[[i]], k=1)
  })
)


sigma_cf_eff0 <- c(
  'tr' = compute_sigma(basis = basis_cf, trt = 1 - data_crossfit$A),
  'nlshrink' = lapply(1:length(degree_of_interactions), function(i) {
    nlshrink_cov((1 - data_crossfit$A) * basis_cf[[i]], k=1)
  })
)




#### HOIF results

## GLM Model
hoif_results <- lapply(names(cf_predictions), function(model_name) {
  preds <- cf_predictions[[model_name]]
  
  eff1 <- compute_HOIF_sequence(
    a_resid   = compute_resid(trt = data_crossfit$A, pred = preds$prob_trt1, type = 'trt'),
    y_resid   = compute_resid(trt = data_crossfit$A, outcome = data_crossfit$Y, pred = preds$prob_outcome1_trt1, type = 'outcome'),
    basis     = rep(basis_cf, length(sigma_cf_eff1)),
    sigma     = sigma_cf_eff1,
    num_cores = 1
  )
  
  eff0 <- compute_HOIF_sequence(
    a_resid   = compute_resid(trt = 1 - data_crossfit$A, pred = 1 - preds$prob_trt1, type = 'trt'),
    y_resid   = compute_resid(trt = 1 - data_crossfit$A, outcome = data_crossfit$Y, pred = preds$prob_outcome1_trt0, type = 'outcome'),
    basis     = rep(basis_cf, length(sigma_cf_eff0)),
    sigma     = sigma_cf_eff0,
    num_cores = 1
  )
  
  list(eff1 = eff1, eff0 = eff0)
})

names(hoif_results) <- names(cf_predictions)


# -------------------------------------------------------
# Results summary: ATE + HOIF corrected ATE per model
# -------------------------------------------------------
results <- lapply(names(cf_predictions), function(model_name) {
  
  aipw <- aipw_results[["glm"]]
  hoif <- hoif_results[["glm"]]
  
  # AIPW ATE
  ate <- aipw$estimate$RD
  
  # HOIF corrected ATE per basis degree
  hoif_corrected <- sapply(1:length(sigma_cf_eff1), function(i) {
    ate - hoif$eff1[2, i] + hoif$eff0[2, i]
  })
  
  # Build summary data frame
  n_rows <- 1 + ncol(hoif$eff1)
  
  data.frame(
    param_num       = rep(param_num, n_rows),
    sim_num         = rep(sim_num, n_rows),
    model           = rep("glm", n_rows),
    estimator       = c("AIPW", paste0("HOIF_basis_dim", hoif$eff1[1, ])),
    basis_dim       = c(NA, hoif$eff1[1, ]),
    sigma_type      = c(NA, colnames(hoif$eff1)),
    AIPW_estimate   = rep(ate, n_rows),
    HOIF_correction = c(NA, hoif$eff1[2, ] - hoif$eff0[2, ]),
    ATE_estimate    = c(ate, hoif_corrected),
    true_ATE        = rep(true_ATE, n_rows),
    bias            = c(ate, hoif_corrected) - true_ATE
  )
})

results_df <- do.call(rbind, results)
rownames(results_df)<-NULL
print(results_df)


output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, 
                         sprintf("results_run%02d_param%02d_sim%03d.RData", run_id, param_num, sim_num))

save(results_df, aipw_results, hoif_results, true_ATE, param_num, sim_num,
     file = output_file)

cat(sprintf("Results saved to: %s\n", output_file))

