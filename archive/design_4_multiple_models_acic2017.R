source('./src/dependencies.R')
source('./src/estimation_functions.R')

# Install/load aciccomp2017 if not already available
if (!requireNamespace("aciccomp2017", quietly = TRUE)) {
  devtools::install_github("vdorie/aciccomp/2017")
}
library(aciccomp2017)

# Design Parameters
param_num <- 1   # aciccomp2017 setting (1:32)
sim_num   <- 1   # aciccomp2017 replication (1:250)

n_oracle  <- 10000
n_train   <- 500
n_est     <- 500

degree_of_interactions <- c(1, 2, 3, 4, 5)

set.seed(123)

##########################################################
## Helper: extract a standardised data.frame from dgp_2017
##   Returns a data.frame with all covariates + A + Y,
##   with column names matching what the rest of the script
##   expects (covariates named as-is from input_2017, plus
##   A for treatment and Y for outcome).
##########################################################
extract_acic_data <- function(param_num, sim_num) {
  df <- dgp_2017(param_num, sim_num)       # returns a data.frame directly
  df <- df %>% rename(A = z, Y = y)
  df
}


# Full covariate names (everything except A and Y)
get_covariate_names <- function(df) {
  names(df)[!names(df) %in% c("A", "Y")]
}


##########################################################
## Oracle data: draw n_oracle rows by repeated calls to
##   dgp_2017 across replications, or simply subsample
##   from the full input_2017 dataset regenerated at scale.
##   Here we use replications 1:n_oracle_reps to build a
##   large pool (each call returns ~4000 rows; we stack and
##   take the first n_oracle rows).
##########################################################

# Each dgp_2017 call returns nrow(input_2017) rows (~4802).
# One call is already > 10 000 / 4802, so we stack two reps.
oracle_raw <- bind_rows(
  extract_acic_data(param_num, 1),
  extract_acic_data(param_num, 2)
)
sim_data_oracle <- oracle_raw[1:n_oracle, ]

oracle_A <- sim_data_oracle$A
oracle_Y <- sim_data_oracle$Y

covariate_names_oracle <- get_covariate_names(sim_data_oracle)

# Compute oracle sigma and basis
basis_oracle <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data       = sim_data_oracle %>% dplyr::select(all_of(covariate_names_oracle)),
    binary_vars = covariate_names_oracle,
    degree_of_interactions = degree_of_interactions[k]
  )
})

sigma_oracle_eff1 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = sim_data_oracle$A))

sigma_oracle_eff0 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = 1 - sim_data_oracle$A))


##########################################################
## Training data: simulation replication sim_num
##   Sample n_train rows from the dgp_2017 output.
##########################################################
dgp_tr_full  <- extract_acic_data(param_num, sim_num)
train_idx    <- sample(nrow(dgp_tr_full), n_train, replace = FALSE)
sim_data_tr  <- dgp_tr_full[train_idx, ]

train_A <- sim_data_tr$A
train_Y <- sim_data_tr$Y

covariate_names_tr <- get_covariate_names(sim_data_tr)

basis_tr <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data        = sim_data_tr %>% dplyr::select(all_of(covariate_names_tr)),
    binary_vars = covariate_names_tr,
    degree_of_interactions = degree_of_interactions[k]
  )
})

sigma_tr_eff1 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = sim_data_tr$A),
                   'nlshrink' = lapply(1:length(degree_of_interactions), function(i) {
                     nlshrink_cov(sim_data_tr$A * basis_tr[[i]], k = 1)
                   }),
                   'unequal_shrink' = lapply(1:length(degree_of_interactions), function(i) {
                     shrinkcovmat.unequal(t(sim_data_tr$A * basis_tr[[i]]), centered = T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(degree_of_interactions), function(i) {
                     shrinkcovmat.equal(t(sim_data_tr$A * basis_tr[[i]]), centered = T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(degree_of_interactions), function(i) {
                     shrinkcovmat.identity(t(sim_data_tr$A * basis_tr[[i]]), centered = T)$Sigmahat
                   }))

sigma_tr_eff0 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = 1 - sim_data_tr$A),
                   'nlshrink' = lapply(1:length(degree_of_interactions), function(i) {
                     nlshrink_cov((1 - sim_data_tr$A) * basis_tr[[i]], k = 1)
                   }),
                   'unequal_shrink' = lapply(1:length(degree_of_interactions), function(i) {
                     shrinkcovmat.unequal(t((1 - sim_data_tr$A) * basis_tr[[i]]), centered = T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(degree_of_interactions), function(i) {
                     shrinkcovmat.equal(t((1 - sim_data_tr$A) * basis_tr[[i]]), centered = T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(degree_of_interactions), function(i) {
                     shrinkcovmat.identity(t((1 - sim_data_tr$A) * basis_tr[[i]]), centered = T)$Sigmahat
                   }))


##########################################################
## Estimation data: use the remaining rows from the same
##   replication (non-overlapping with training).
##########################################################
est_idx      <- setdiff(seq_len(nrow(dgp_tr_full)), train_idx)
est_idx      <- sample(est_idx, min(n_est, length(est_idx)), replace = FALSE)
sim_data_est <- dgp_tr_full[est_idx, ]

covariate_names_est <- get_covariate_names(sim_data_est)

basis_est <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data        = sim_data_est %>% dplyr::select(all_of(covariate_names_est)),
    binary_vars = covariate_names_est,
    degree_of_interactions = degree_of_interactions[k]
  )
})


### *******************************************
### Creating basis from decision tree model
### *******************************************
library(rpart)
library(partykit)

tree_basis_model <- rpart(
  A ~ .,
  data    = sim_data_tr %>% dplyr::select(-Y),
  method  = "class",
  control = rpart.control(
    maxdepth = 4,
    cp       = 0.001,
    minsplit = 20
  )
)

create_tree_basis <- function(tree_model, train_covariates, new_covariates) {
  party_model   <- as.party(tree_model)
  leaf_train    <- as.integer(predict(party_model, newdata = train_covariates, type = "node"))
  terminal_nodes <- sort(unique(leaf_train))
  leaf_new      <- as.integer(predict(party_model, newdata = new_covariates, type = "node"))
  basis         <- outer(leaf_new, terminal_nodes, FUN = "==") * 1
  colnames(basis) <- paste0("leaf_", terminal_nodes)
  as.matrix(basis)
}

X_tr  <- sim_data_tr  %>% dplyr::select(-A, -Y)
X_est <- sim_data_est %>% dplyr::select(-A, -Y)

basis_tree_tr  <- create_tree_basis(tree_basis_model, X_tr, X_tr)
basis_tree_est <- create_tree_basis(tree_basis_model, X_tr, X_est)

basis_tree_est_list    <- list()
basis_tree_est_list[[1]] <- basis_tree_est

# Sigma on estimation data
sigma_est_eff1 <- c('est' = compute_sigma(basis = basis_est,
                                          trt = sim_data_est$A))
sigma_est_eff0 <- c('est' = compute_sigma(basis = basis_est,
                                          trt = 1 - sim_data_est$A))

sigma_est_tree_eff1 <- c('est' = compute_sigma(basis = basis_tree_est_list,
                                               trt = sim_data_est$A))
sigma_est_tree_eff0 <- c('est' = compute_sigma(basis = basis_tree_est_list,
                                               trt = 1 - sim_data_est$A))


##########################################################
## Estimate nuisance parameter models using training data
##########################################################

# GLM Model
params_glm_trt     <- find_params_glm(binary_vars = covariate_names_tr)
params_glm_outcome <- find_params_glm(binary_vars = c(covariate_names_tr, "A"))

# KNN Model
params_knn_trt <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                  label_vector  = sim_data_tr$A,
                                  nfold = 5,
                                  k     = seq(11, 101, 2),
                                  num_cores = 10)
params_knn_outcome <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(-Y),
                                      label_vector  = sim_data_tr$Y,
                                      nfold = 5,
                                      k     = seq(11, 101, 2),
                                      num_cores = 10)

# Lasso Model
params_lasso_trt <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                      label_vector  = sim_data_tr$A,
                                      binary_vars   = covariate_names_tr,
                                      degree_of_interactions = 3,
                                      nfold = 5)
params_lasso_outcome <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(-Y),
                                          label_vector  = sim_data_tr$Y,
                                          binary_vars   = c(covariate_names_tr, "A"),
                                          degree_of_interactions = 3,
                                          nfold = 5)

# Random Forest Model
params_random_forest_trt <- find_params_random_forest_model(
  covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
  label_vector  = sim_data_tr$A,
  nfold     = 5,
  num_trees = seq(500, 1000, 50),
  num_vars  = seq(1, 5, 1),
  num_cores = 10)
params_random_forest_outcome <- find_params_random_forest_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  nfold     = 5,
  num_trees = seq(500, 1000, 50),
  num_vars  = seq(1, 5, 1),
  num_cores = 10)


##########################################################
## Predict nuisance parameters on estimation data
##########################################################

## GLM Model
prob_trt1_glm <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
  label_vector  = sim_data_tr$A,
  params        = params_glm_trt,
  model         = "glm",
  predict_data  = sim_data_est %>% dplyr::select(-c(Y, A)))

prob_outcome1_trt1_glm <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_glm_outcome,
  model         = "glm",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 1))

prob_outcome1_trt0_glm <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_glm_outcome,
  model         = "glm",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 0))

## KNN Model
prob_trt1_knn <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
  label_vector  = sim_data_tr$A,
  params        = params_knn_trt,
  model         = "knn",
  predict_data  = sim_data_est %>% dplyr::select(-c(Y, A)))

prob_outcome1_trt1_knn <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_knn_outcome,
  model         = "knn",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 1))

prob_outcome1_trt0_knn <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_knn_outcome,
  model         = "knn",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 0))

## Lasso Model
prob_trt1_lasso <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
  label_vector  = sim_data_tr$A,
  params        = params_lasso_trt,
  model         = "lasso",
  predict_data  = sim_data_est %>% dplyr::select(-c(Y, A)))

prob_outcome1_trt1_lasso <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_lasso_outcome,
  model         = "lasso",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 1))

prob_outcome1_trt0_lasso <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_lasso_outcome,
  model         = "lasso",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 0))

## Random Forest Model
prob_trt1_random_forest <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
  label_vector  = sim_data_tr$A,
  params        = params_random_forest_trt,
  model         = "random_forest",
  predict_data  = sim_data_est %>% dplyr::select(-c(Y, A)))

prob_outcome1_trt1_random_forest <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_random_forest_outcome,
  model         = "random_forest",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 1))

prob_outcome1_trt0_random_forest <- estimate_prob_individual_model(
  covariates_df = sim_data_tr %>% dplyr::select(-Y),
  label_vector  = sim_data_tr$Y,
  params        = params_random_forest_outcome,
  model         = "random_forest",
  predict_data  = sim_data_est %>% dplyr::select(-Y) %>% mutate(A = 0))


#######################
### First Order Estimator (AIPW)
########################

## GLM
estimate_AIPW_glm <- estimate_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_glm,
  prob_outcome1_trt1   = prob_outcome1_trt1_glm,
  prob_outcome1_trt0   = prob_outcome1_trt0_glm)

estimate_var_AIPW_glm <- estimate_var_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_glm,
  prob_outcome1_trt1   = prob_outcome1_trt1_glm,
  prob_outcome1_trt0   = prob_outcome1_trt0_glm)

## KNN
estimate_AIPW_knn <- estimate_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_knn,
  prob_outcome1_trt1   = prob_outcome1_trt1_knn,
  prob_outcome1_trt0   = prob_outcome1_trt0_knn)

estimate_var_AIPW_knn <- estimate_var_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_knn,
  prob_outcome1_trt1   = prob_outcome1_trt1_knn,
  prob_outcome1_trt0   = prob_outcome1_trt0_knn)

## Lasso
estimate_AIPW_lasso <- estimate_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_lasso,
  prob_outcome1_trt1   = prob_outcome1_trt1_lasso,
  prob_outcome1_trt0   = prob_outcome1_trt0_lasso)

estimate_var_AIPW_lasso <- estimate_var_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_lasso,
  prob_outcome1_trt1   = prob_outcome1_trt1_lasso,
  prob_outcome1_trt0   = prob_outcome1_trt0_lasso)

## Random Forest
estimate_AIPW_random_forest <- estimate_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_random_forest,
  prob_outcome1_trt1   = prob_outcome1_trt1_random_forest,
  prob_outcome1_trt0   = prob_outcome1_trt0_random_forest)

estimate_var_AIPW_random_forest <- estimate_var_AIPW(
  trt                  = sim_data_est$A,
  outcome              = sim_data_est$Y,
  prob_trt1            = prob_trt1_random_forest,
  prob_outcome1_trt1   = prob_outcome1_trt1_random_forest,
  prob_outcome1_trt0   = prob_outcome1_trt0_random_forest)


#######################
### Second Order Estimator (HOIF)
########################

## GLM
HOIF_glm_eff1 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = sim_data_est$A, pred = prob_trt1_glm, type = 'trt'),
  y_resid  = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_glm, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1, sigma_est_tree_eff1),
  num_cores = 1)

HOIF_glm_eff0 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = 1 - sim_data_est$A, pred = 1 - prob_trt1_glm, type = 'trt'),
  y_resid  = compute_resid(trt = 1 - sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_glm, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1)

## KNN
HOIF_knn_eff1 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = sim_data_est$A, pred = prob_trt1_knn, type = 'trt'),
  y_resid  = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_knn, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1, sigma_est_tree_eff1),  # fixed: eff1
  num_cores = 1)

HOIF_knn_eff0 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = 1 - sim_data_est$A, pred = 1 - prob_trt1_knn, type = 'trt'),
  y_resid  = compute_resid(trt = 1 - sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_knn, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1)

## Lasso
HOIF_lasso_eff1 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = sim_data_est$A, pred = prob_trt1_lasso, type = 'trt'),
  y_resid  = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_lasso, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1, sigma_est_tree_eff1),  # fixed: eff1
  num_cores = 1)

HOIF_lasso_eff0 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = 1 - sim_data_est$A, pred = 1 - prob_trt1_lasso, type = 'trt'),
  y_resid  = compute_resid(trt = 1 - sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_lasso, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1)

## Random Forest
HOIF_random_forest_eff1 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = sim_data_est$A, pred = prob_trt1_random_forest, type = 'trt'),
  y_resid  = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_random_forest, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1, sigma_est_tree_eff1),  # fixed: eff1
  num_cores = 1)

HOIF_random_forest_eff0 <- compute_HOIF_sequence(
  a_resid  = compute_resid(trt = 1 - sim_data_est$A, pred = 1 - prob_trt1_random_forest, type = 'trt'),
  y_resid  = compute_resid(trt = 1 - sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_random_forest, type = 'outcome'),
  basis    = rep(basis_est, 8),
  sigma    = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1)
