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


### ML kernel

### ML kernel code 

d_idx <- 2   # corresponds to degree_of_interactions[2]

X_tr <- sim_data_tr %>% dplyr::select(-A, -Y) %>% as.matrix()

Phi_tr <- as.matrix(basis_tr[[d_idx]])

# choose which covariance estimate you want
Sigma_eff1 <- sigma_tr_eff1[[paste0("tr", d_idx)]]
Sigma_eff0 <- sigma_tr_eff0[[paste0("tr", d_idx)]]

# regularized inverse for stability
ridge <- 1e-6

Sigma_eff1_inv <- solve(Sigma_eff1 + ridge * diag(ncol(Sigma_eff1)))
Sigma_eff0_inv <- solve(Sigma_eff0 + ridge * diag(ncol(Sigma_eff0)))

# training kernel matrices for HOIF
K_eff1 <- Phi_tr %*% Sigma_eff1_inv %*% t(Phi_tr)
K_eff0 <- Phi_tr %*% Sigma_eff0_inv %*% t(Phi_tr)

### Creating the train data
X_tr <- sim_data_tr %>% dplyr::select(-A, -Y) %>% as.matrix()

make_pair_features <- function(X1, X2) {
  out <- cbind(X1, X2, abs(X1 - X2), X1 * X2)
  
  p <- ncol(X1)
  colnames(out) <- c(
    paste0("x1_", colnames(X1)),
    paste0("x2_", colnames(X2)),
    paste0("absdiff_", colnames(X1)),
    paste0("prod_", colnames(X1))
  )
  
  out
}

pair_idx <- which(upper.tri(K_eff1), arr.ind = TRUE)

set.seed(123)
sel <- sample(1:nrow(pair_idx), size = min(20000, nrow(pair_idx)))

idx_i <- pair_idx[sel, 1]
idx_j <- pair_idx[sel, 2]

pair_X <- make_pair_features(X_tr[idx_i, ], X_tr[idx_j, ])
pair_y <- K_eff1[cbind(idx_i, idx_j)]


### Fit a decision tree model
pair_df <- as.data.frame(pair_X)
pair_df$y <- pair_y

set.seed(123)
train_ids <- sample(seq_len(nrow(pair_df)), size = floor(0.8 * nrow(pair_df)))

train_df <- pair_df[train_ids, ]
test_df  <- pair_df[-train_ids, ]

tree_model <- rpart(
  y ~ .,
  data = train_df,
  method = "anova",
  control = rpart.control(
    maxdepth = ,
    cp = 0.001,
    minsplit = 20
  )
)

pred <- predict(tree_model, newdata = test_df)

mean((pred - test_df$y)^2)
cor(pred, test_df$y)
sqrt(mean((pred - test_df$y)^2)) / sd(test_df$y)

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

sigma_est_eff1 <- c('est' = compute_sigma(basis = basis_est,
                                          trt = sim_data_est$A))

sigma_est_eff0 <- c('est' = compute_sigma(basis = basis_est,
                                          trt = 1-sim_data_est$A))


### ML estimation
d_idx <- 2

Phi_est <- as.matrix(basis_est[[d_idx]])
Sigma_est_eff1 <- sigma_est_eff1[[paste0("est", d_idx)]]
Sigma_est_eff0 <- sigma_est_eff0[[paste0("est", d_idx)]]
ridge <- 1e-6
Sigma_est_eff1_inv <- solve(Sigma_est_eff1 + ridge * diag(ncol(Sigma_est_eff1)))
Sigma_est_eff0_inv <- solve(Sigma_est_eff0 + ridge * diag(ncol(Sigma_est_eff0)))
K_est_eff1 <- Phi_est %*% Sigma_est_eff1_inv %*% t(Phi_est)
K_est_eff0 <- Phi_est %*% Sigma_est_eff0_inv %*% t(Phi_est)


## ML estimatioin
X_est <- sim_data_est %>% dplyr::select(-A, -Y) %>% as.matrix()
n_est <- nrow(X_est)

Khat_est <- matrix(0, n_est, n_est)

for (i in 1:n_est) {
  X1 <- X_est[rep(i, n_est), , drop = FALSE]
  X2 <- X_est
  
  pair_feat <- make_pair_features(X1, X2)
  pair_feat <- as.data.frame(pair_feat)
  
  Khat_est[i, ] <- predict(tree_model, newdata = pair_feat)
}

# enforce symmetry
Khat_est <- 0.5 * (Khat_est + t(Khat_est))

cor(as.vector(Khat_est), as.vector(K_est_eff1))
mean((Khat_est - K_est_eff1)^2)
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

compute_HOIF_22_from_kernel <- function(y_resid, a_resid, Khat) {
  
  n <- length(y_resid)
  
  stat <- ((as.numeric(t(y_resid) %*% Khat %*% a_resid) -
              sum(y_resid * diag(Khat) * a_resid)) / (n * (n - 1))) * (-1)
  
  return(stat)
}

## GLM Model
HOIF_glm_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_glm, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_glm, type = 'outcome'),
  basis = rep(basis_est,7),
  sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
  num_cores = 1
)

HOIF_glm_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_glm, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_glm, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
  num_cores = 1
)


## KNN Model
HOIF_knn_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_knn, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_knn, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
  num_cores = 1
)

HOIF_knn_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_knn, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_knn, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
  num_cores = 1
)


## Lasso Model
HOIF_lasso_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_lasso, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_lasso, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
  num_cores = 1
)

HOIF_lasso_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_lasso, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_lasso, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
  num_cores = 1
)


## Random Forest Model
HOIF_random_forest_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_random_forest, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_random_forest, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff1, sigma_tr_eff1, sigma_est_eff1),
  num_cores = 1
)

HOIF_random_forest_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_random_forest, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_random_forest, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
  num_cores = 1
)

HOIF_random_forest_kernel <- compute_HOIF_22_from_kernel(
  a_resid = compute_resid(trt = 1-sim_data_est$A, pred = 1-prob_trt1_random_forest, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_random_forest, type = 'outcome'),
  Khat = Khat_est
)
