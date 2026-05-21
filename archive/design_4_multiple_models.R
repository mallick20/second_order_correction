source('./src/dependencies.R')
source('./src/estimation_functions.R')

# Design Parameters
n_oracle <- 10000
n_train <- 500
n_est <- 500

## For continuous covariates
sigma2 <- 1
degree_of_interactions <- c(1,2,3,4,5)

set.seed(123)

# Construct the treatment dependent and outcome based only on covariates - no treatment effect
key_trt <- expand.grid(Z1=c(0,1), Z2=c(0,1), Z3=c(0,1), Z4=c(0,1)) %>%
  {bind_cols(., val_trt=sample(0.05*1:nrow(.), replace=F))}

key_outcome <- expand.grid(Z1=c(0,1), Z2=c(0,1), Z3=c(0,1), Z4=c(0,1)) %>%
  {bind_cols(., val_outcome=sample(0.05*1:nrow(.), replace=F))}


##########################################################
## simulate oracle data and parameters
##########################################################
sim_data_oracle <- data.frame(
  Z1 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z2 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z3 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z4 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z5 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z6 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z7 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z8 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z9 = rnorm(n_oracle, mean = 0, sd = sqrt(sigma2)),
  Z10 = rnorm(n_oracle, mean = 0, sd = sqrt(sigma2))
) %>% 
  merge(., key_trt) %>%
  merge(., key_outcome) %>% 
  mutate(
    A = rbinom(n(), 1, val_trt),
    Y = rbinom(n(), 1, val_outcome)
  ) %>% dplyr::select(-c(val_trt, val_outcome))

oracle_A <- sim_data_oracle$A
oracle_Y <- sim_data_oracle$Y

# Compute oracle sigma and basis
basis_oracle <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data = sim_data_oracle %>% dplyr::select(-c(A,Y)),
    binary_vars = names(sim_data_oracle %>% dplyr::select(-c(A,Y))),
    degree_of_interactions = degree_of_interactions[k]
  )
})


sigma_oracle_eff1 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = sim_data_oracle$A))

sigma_oracle_eff0 <- c('oracle' = compute_sigma(basis = basis_oracle,
                                                trt = 1-sim_data_oracle$A))




##########################################################
## simulate training data
##########################################################
sim_data_tr <- data.frame(
  Z1 = as.numeric(rbernoulli(n_train,p=0.5)), Z2 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z3 = as.numeric(rbernoulli(n_train,p=0.5)), Z4 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z5 = as.numeric(rbernoulli(n_train,p=0.5)), Z6 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z7 = as.numeric(rbernoulli(n_train,p=0.5)), Z8 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z9 = rnorm(n_train, mean = 0, sd = sqrt(sigma2)),
  Z10 = rnorm(n_train, mean = 0, sd = sqrt(sigma2))
) %>%
  merge(., key_trt) %>%
  merge(., key_outcome) %>% 
  mutate(
    A = rbinom(n(), 1, val_trt),      # Treatment based on Z1-z4
    Y = rbinom(n(), 1, val_outcome)   # Outcome based on Z1-z4
  ) 

train_A <- sim_data_tr$A
train_Y <- sim_data_tr$Y

sim_data_tr <- sim_data_tr %>% dplyr::select(-c(val_trt,val_outcome))

basis_tr <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data = sim_data_tr %>% dplyr::select(-c(A,Y)),
    binary_vars = names(sim_data_tr %>% dplyr::select(-c(A,Y))),
    degree_of_interactions = degree_of_interactions[k]
  )
})

sigma_tr_eff1 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = sim_data_tr$A),
                   'nlshrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     nlshrink_cov(sim_data_tr$A*basis_tr[[i]], k=1)
                   }),
                   'unequal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.unequal(t(sim_data_tr$A*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.equal(t(sim_data_tr$A*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.identity(t(sim_data_tr$A*basis_tr[[i]]), centered=T)$Sigmahat
                   }))

sigma_tr_eff0 <- c('tr' = compute_sigma(basis = basis_tr,
                                        trt = 1-sim_data_tr$A),
                   'nlshrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     nlshrink_cov((1-sim_data_tr$A)*basis_tr[[i]], k=1)
                   }),
                   'unequal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.unequal(t((1-sim_data_tr$A)*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'equal_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.equal(t((1-sim_data_tr$A)*basis_tr[[i]]), centered=T)$Sigmahat
                   }),
                   'identity_shrink' = lapply(1:length(c(degree_of_interactions)), function(i) {
                     shrinkcovmat.identity(t((1-sim_data_tr$A)*basis_tr[[i]]), centered=T)$Sigmahat
                   }))


### Estimated data
n = n_est
sim_data_est <- data.frame(
  Z1 = as.numeric(rbernoulli(n,p=0.5)), Z2 = as.numeric(rbernoulli(n,p=0.5)),
  Z3 = as.numeric(rbernoulli(n,p=0.5)), Z4 = as.numeric(rbernoulli(n,p=0.5)),
  Z5 = as.numeric(rbernoulli(n,p=0.5)), Z6 = as.numeric(rbernoulli(n,p=0.5)),
  Z7 = as.numeric(rbernoulli(n,p=0.5)), Z8 = as.numeric(rbernoulli(n,p=0.5)),
  Z9 = rnorm(n, mean = 0, sd = sqrt(sigma2)),
  Z10 = rnorm(n, mean = 0, sd = sqrt(sigma2))
) %>% 
  merge(., key_trt) %>%
  merge(., key_outcome) %>% 
  mutate(
    A = rbinom(n(), 1, val_trt),
    Y = rbinom(n(), 1, val_outcome)
  ) 

p_A <- sim_data_est$val_trt
p_Y <- sim_data_est$val_outcome

sim_data_est <- sim_data_est %>% dplyr::select(-c(val_trt,val_outcome))

basis_est <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data = sim_data_est %>% dplyr::select(-c(A,Y)),
    binary_vars = names(sim_data_est %>% dplyr::select(-c(A,Y))),
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
  data = sim_data_tr %>% dplyr::select(-Y),
  method = "class",
  control = rpart.control(
    maxdepth = 4,
    cp = 0.001,
    minsplit = 20
  )
)



create_tree_basis <- function(tree_model, train_covariates, new_covariates) {
  
  party_model <- as.party(tree_model)
  
  ## leaf ids on train data define the basis columns
  leaf_train <- predict(party_model, newdata = train_covariates, type = "node")
  leaf_train <- as.integer(leaf_train)
  terminal_nodes <- sort(unique(leaf_train))
  
  ## leaf ids on new data
  leaf_new <- predict(party_model, newdata = new_covariates, type = "node")
  leaf_new <- as.integer(leaf_new)
  
  ## one-hot basis
  basis <- outer(leaf_new, terminal_nodes, FUN = "==") * 1
  
  colnames(basis) <- paste0("leaf_", terminal_nodes)
  basis <- as.matrix(basis)
  
  return(basis)
}


X_tr  <- sim_data_tr  %>% dplyr::select(-A, -Y)
X_est <- sim_data_est %>% dplyr::select(-A, -Y)

basis_tree_tr <- create_tree_basis(
  tree_model = tree_basis_model,
  train_covariates = X_tr,
  new_covariates = X_tr
)

basis_tree_est <- create_tree_basis(
  tree_model = tree_basis_model,
  train_covariates = X_tr,
  new_covariates = X_est
)

basis_tree_est_list <- list()
basis_tree_est_list[[1]] <- basis_tree_est


# Calculating sigma
sigma_est_eff1 <- c('est' = compute_sigma(basis = basis_est,
                                          trt = sim_data_est$A))

sigma_est_eff0 <- c('est' = compute_sigma(basis = basis_est,
                                          trt = 1-sim_data_est$A))



sigma_est_tree_eff1 <- c('est' = compute_sigma(basis = basis_tree_est_list,
                                          trt = sim_data_est$A))

sigma_est_tree_eff0 <- c('est' = compute_sigma(basis = basis_tree_est_list,
                                          trt = 1-sim_data_est$A))


##########################################################
## estimate nuisance parameter models using training data
##########################################################

# GLM Model
params_glm_trt <- find_params_glm(binary_vars = names(sim_data_tr %>% dplyr::select(-c(A,Y))))
params_glm_outcome <- find_params_glm(binary_vars = names(sim_data_tr %>% dplyr::select(-c(Y))))

# KNN Model
params_knn_trt <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                  label_vector = sim_data_tr %>% dplyr::select(A) %>% {.[[1]]},
                                  nfold = 5,
                                  k = seq(11,101,2),
                                  num_cores = 10)
params_knn_outcome <- find_params_knn(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                      label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                      nfold = 5,
                                      k = seq(11,101,2),
                                      num_cores = 10)

# Lasso Model
params_lasso_trt <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                      label_vector = sim_data_tr %>% dplyr::select(A) %>% {.[[1]]},
                                      binary_vars = names(sim_data_tr %>% dplyr::select(-c(A,Y))),
                                      degree_of_interactions = 3,
                                      nfold = 5)
params_lasso_outcome <- find_params_lasso(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                          label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                          binary_vars = names(sim_data_tr %>% dplyr::select(-c(Y))),
                                          degree_of_interactions = 3, 
                                          nfold = 5)

# Random Forest Model
params_random_forest_trt <- find_params_random_forest_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                                            label_vector = sim_data_tr %>% dplyr::select(A) %>% {.[[1]]},
                                                            nfold = 5,
                                                            num_trees = seq(500,1000,50),
                                                            num_vars = seq(1,5,1),
                                                            num_cores = 10)
params_random_forest_outcome <- find_params_random_forest_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                nfold = 5,
                                                                num_trees = seq(500,1000,50),
                                                                num_vars = seq(1,5,1),
                                                                num_cores = 10)




##########################################################
## predict nuisance parameters using estimation data
##########################################################

## GLM Model
# Predict treatment (Propensity)
prob_trt1_glm <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                                label_vector = sim_data_tr %>% dplyr::select(A) %>% {.[[1]]},
                                                params = params_glm_trt, 
                                                model = "glm",
                                                predict_data = sim_data_est %>% dplyr::select(-c(Y, A)))

# Outcome regression models (Outcome probability when treatment and control)
prob_outcome1_trt1_glm <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                         params = params_glm_outcome, 
                                                         model = "glm",
                                                         predict_data = sim_data_tr %>% dplyr::select(-c(Y)) %>% mutate(A = 1))
prob_outcome1_trt0_glm <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                         params = params_glm_outcome, 
                                                         model = "glm",
                                                         predict_data = sim_data_tr %>% dplyr::select(-c(Y)) %>% mutate(A = 0))


## KNN Model
prob_trt1_knn <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                                label_vector = sim_data_tr %>% dplyr::select(A) %>% {.[[1]]},
                                                params = params_knn_trt, 
                                                model = "knn",
                                                predict_data = sim_data_est %>% dplyr::select(-c(Y, A)))

prob_outcome1_trt1_knn <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                         params = params_knn_outcome, 
                                                         model = "knn",
                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(A = 1))

prob_outcome1_trt0_knn <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                         label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                         params = params_knn_outcome, 
                                                         model = "knn",
                                                         predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(A = 0))


## Lasso Model
prob_trt1_lasso <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                                  label_vector = sim_data_tr %>% dplyr::select(A) %>% {.[[1]]},
                                                  params = params_lasso_trt, 
                                                  model = "lasso",
                                                  predict_data = sim_data_est %>% dplyr::select(-c(Y, A)))

prob_outcome1_trt1_lasso <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                           label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                           params = params_lasso_outcome, 
                                                           model = "lasso",
                                                           predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(A = 1))

prob_outcome1_trt0_lasso <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                           label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                           params = params_lasso_outcome, 
                                                           model = "lasso",
                                                           predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(A = 0))


## Random Forest Model
prob_trt1_random_forest <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y, A)),
                                                          label_vector = sim_data_tr %>% dplyr::select(A) %>% {.[[1]]},
                                                          params = params_random_forest_trt, 
                                                          model = "random_forest",
                                                          predict_data = sim_data_est %>% dplyr::select(-c(Y, A)))

prob_outcome1_trt1_random_forest <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                   label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                   params = params_random_forest_outcome, 
                                                                   model = "random_forest",
                                                                   predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(A = 1))

prob_outcome1_trt0_random_forest <- estimate_prob_individual_model(covariates_df = sim_data_tr %>% dplyr::select(-c(Y)),
                                                                   label_vector = sim_data_tr %>% dplyr::select(Y) %>% {.[[1]]},
                                                                   params = params_random_forest_outcome, 
                                                                   model = "random_forest",
                                                                   predict_data = sim_data_est %>% dplyr::select(-c(Y)) %>% mutate(A = 0))

#######################
### First Order Estimator
########################

## GLM Model
estimate_AIPW_glm <- 
  estimate_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_glm,
    prob_outcome1_trt1 = prob_outcome1_trt1_glm,
    prob_outcome1_trt0 = prob_outcome1_trt0_glm
  )

estimate_var_AIPW_glm <- 
  estimate_var_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_glm,
    prob_outcome1_trt1 = prob_outcome1_trt1_glm,
    prob_outcome1_trt0 = prob_outcome1_trt0_glm
  )

## KNN Model
estimate_AIPW_knn <- 
  estimate_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_knn,
    prob_outcome1_trt1 = prob_outcome1_trt1_knn,
    prob_outcome1_trt0 = prob_outcome1_trt0_knn
  )

estimate_var_AIPW_knn <- 
  estimate_var_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_knn,
    prob_outcome1_trt1 = prob_outcome1_trt1_knn,
    prob_outcome1_trt0 = prob_outcome1_trt0_knn
  )

## Lasso Model
estimate_AIPW_lasso <- 
  estimate_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_lasso,
    prob_outcome1_trt1 = prob_outcome1_trt1_lasso,
    prob_outcome1_trt0 = prob_outcome1_trt0_lasso
  )

estimate_var_AIPW_lasso <- 
  estimate_var_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_lasso,
    prob_outcome1_trt1 = prob_outcome1_trt1_lasso,
    prob_outcome1_trt0 = prob_outcome1_trt0_lasso
  )

## Random Forest Model
estimate_AIPW_random_forest <- 
  estimate_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_random_forest,
    prob_outcome1_trt1 = prob_outcome1_trt1_random_forest,
    prob_outcome1_trt0 = prob_outcome1_trt0_random_forest
  )

estimate_var_AIPW_random_forest <- 
  estimate_var_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_random_forest,
    prob_outcome1_trt1 = prob_outcome1_trt1_random_forest,
    prob_outcome1_trt0 = prob_outcome1_trt0_random_forest
  )

#######################
### Second order estimator
########################

## GLM Model
HOIF_glm_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_glm, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_glm, type = 'outcome'),
  basis = rep(basis_est,8),
  sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1, sigma_est_tree_eff1),
  num_cores = 1
)

HOIF_glm_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_glm, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_glm, type = 'outcome'),
  basis = rep(basis_est, 8),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1
)


## KNN Model
HOIF_knn_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_knn, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_knn, type = 'outcome'),
  basis = rep(basis_est, 8),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1
)

HOIF_knn_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_knn, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_knn, type = 'outcome'),
  basis = rep(basis_est, 8),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1
)


## Lasso Model
HOIF_lasso_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_lasso, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_lasso, type = 'outcome'),
  basis = rep(basis_est, 8),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1
)

HOIF_lasso_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_lasso, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_lasso, type = 'outcome'),
  basis = rep(basis_est, 8),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1
)


## Random Forest Model
HOIF_random_forest_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_random_forest, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_random_forest, type = 'outcome'),
  basis = rep(basis_est, 8),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1
)

HOIF_random_forest_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_random_forest, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_random_forest, type = 'outcome'),
  basis = rep(basis_est, 8),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0, sigma_est_tree_eff0),
  num_cores = 1
)

