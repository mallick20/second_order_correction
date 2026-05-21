source('./src/dependencies.R')
source('./src/estimation_functions.R')

# Basic design which consist of few parameters and a continuous treatment variable.
n_train <- 500
n_oracle <- 100000
degree_of_interactions <- c(1,2,3,4,5,6,7)

set.seed(123)

# Construct the treatment dependent on only one covariate and another is a noise variable
# Random treatment
# parameters for true function
key_trt <- expand.grid(Z1=c(0,1), Z2=c(0,1), Z3=c(0,1), Z4=c(0,1)) %>%
  {bind_cols(., val_trt=sample(0.05*1:nrow(.), replace=F))}

key_outcome <- expand.grid(Z1=c(0,1), Z2=c(0,1), Z3=c(0,1), Z4=c(0,1)) %>%
  {bind_cols(., val_outcome=sample(0.05*1:nrow(.), replace=F))}

# There is this one confounder and we will create only one more noise variable.
# This is the oracle data
sim_data_oracle <- data.frame(
  Z1 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z2 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z3 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z4 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z5 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z6 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z7 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z8 = as.numeric(rbernoulli(n_oracle,p=0.5)),
  Z9 = as.numeric(rbernoulli(n_oracle,p=0.5)), Z10 = as.numeric(rbernoulli(n_oracle,p=0.5))
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


## Simple ATE calculation 
## 1. expected value of y given treatment - control
ate_oracle <- 0

##########################################################
## simulate training data
##########################################################
sim_data_train <- data.frame(
  Z1 = as.numeric(rbernoulli(n_train,p=0.5)), Z2 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z3 = as.numeric(rbernoulli(n_train,p=0.5)), Z4 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z5 = as.numeric(rbernoulli(n_train,p=0.5)), Z6 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z7 = as.numeric(rbernoulli(n_train,p=0.5)), Z8 = as.numeric(rbernoulli(n_train,p=0.5)),
  Z9 = as.numeric(rbernoulli(n_train,p=0.5)), Z10 = as.numeric(rbernoulli(n_train,p=0.5))
) %>%
  merge(., key_trt) %>%
  merge(., key_outcome) %>% 
  mutate(
    A = rbinom(n(), 1, val_trt),      # Treatment based on Z1-z4
    Y = rbinom(n(), 1, val_outcome)   # Outcome based on Z1-z4
  ) 

train_A <- sim_data_train$A
train_Y <- sim_data_train$Y

sim_data_train <- sim_data_train %>% dplyr::select(-c(val_trt,val_outcome))

basis_tr <- lapply(1:length(degree_of_interactions), function(k) {
  create_binary_var_basis(
    data = sim_data_train %>% dplyr::select(-c(A,Y)),
    binary_vars = names(sim_data_train %>% dplyr::select(-c(A,Y))),
    degree_of_interactions = degree_of_interactions[k]
  )
})

sigma_train_eff1 <- c('tr' = compute_sigma(basis = basis_tr,
                                           trt = sim_data_train$A))

sigma_train_eff0 <- c('tr' = compute_sigma(basis = basis_tr,
                                           trt = 1-sim_data_train$A))

### Estimated data
n = 500
sim_data_est <- data.frame(
  Z1 = as.numeric(rbernoulli(n,p=0.5)), Z2 = as.numeric(rbernoulli(n,p=0.5)),
  Z3 = as.numeric(rbernoulli(n,p=0.5)), Z4 = as.numeric(rbernoulli(n,p=0.5)),
  Z5 = as.numeric(rbernoulli(n,p=0.5)), Z6 = as.numeric(rbernoulli(n,p=0.5)),
  Z7 = as.numeric(rbernoulli(n,p=0.5)), Z8 = as.numeric(rbernoulli(n,p=0.5)),
  Z9 = as.numeric(rbernoulli(n,p=0.5)), Z10 = as.numeric(rbernoulli(n,p=0.5))
) %>% 
  merge(., key_trt) %>%
  merge(., key_outcome) %>% 
  mutate(
    A = rbinom(n(), 1, val_trt),
    Y = rbinom(n(), 1, val_outcome)
  ) 

p_X <- sim_data_est$val_trt
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

### AIPW ATE
# Doubly robust (AIPW) estimator for ATE: E[Y(1) - Y(0)]
# Assumes unconfoundedness given Z1 (and optionally Z2)

aipw_ate <- function(dat, outcome = "Y", treat = "A", covars = c("Z1")) {
  Y <- dat[[outcome]]
  A <- dat[[treat]]
  
  # 1) Propensity score model: e(x) = P(A=1 | X)
  f_e <- as.formula(paste("A ~", paste(covars, collapse = " + ")))
  m_e <- glm(f_e, data = dat, family = binomial())
  ehat <- as.numeric(predict(m_e, type = "response"))
  
  # guard against 0/1 propensity
  eps <- 1e-6
  ehat <- pmin(pmax(ehat, eps), 1 - eps)
  
  # 2) Outcome regression: mu(a,x) = E[Y | A=a, X]
  f_mu <- as.formula(paste("Y ~ A ", "+", paste(covars, collapse = " + ")))
  m_mu <- glm(f_mu, data = dat, family = binomial())
  
  # Predict mu1(x) and mu0(x)
  dat1 <- dat; dat1[[treat]] <- 1
  dat0 <- dat; dat0[[treat]] <- 0
  mu1 <- as.numeric(predict(m_mu, newdata = dat1, type = "response"))
  mu0 <- as.numeric(predict(m_mu, newdata = dat0, type = "response"))
  
  # 3) AIPW / Doubly-robust score and ATE
  psi <- (mu1 - mu0) + (A * (Y - mu1) / ehat) - ((1 - A) * (Y - mu0) / (1 - ehat))
  ate_hat <- mean(psi)
  
  # Simple SE via influence-function (large-sample)
  se_hat <- sd(psi) / sqrt(nrow(dat))
  
  list(ate = ate_hat, se = se_hat, e_model = m_e, mu_model = m_mu)
}

# --- Use it on your simulated training data ---
res <- aipw_ate(sim_data_train, outcome = "Y", treat = "A", covars = c("Z1", "Z2", "Z3", "Z4"))
ate_train_aipw <- res$ate
res$se
c(res$ate - 1.96 * res$se, res$ate + 1.96 * res$se)

res_2 <- aipw_ate(sim_data_train, outcome = "Y", treat = "A", covars = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8", "Z9", "Z10"))
ate_train_aipw_2 <- res_2$ate
res_2$se
c(res_2$ate - 1.96 * res_2$se, res_2$ate + 1.96 * res_2$se)

#######################
### Second order estimator
########################




#### Just use glm model for now
params_glm_trt <- find_params_glm(binary_vars = names(sim_data_train %>% dplyr::select(-c(A,Y))))
params_glm_outcome <- find_params_glm(binary_vars = names(sim_data_train %>% dplyr::select(-c(Y))))

# Predict treatment (Propensity)
prob_trt1_glm <- estimate_prob_individual_model(covariates_df = sim_data_train %>% dplyr::select(-c(Y, A)),
                                                label_vector = sim_data_train %>% dplyr::select(A) %>% {.[[1]]},
                                                params = params_glm_trt, 
                                                model = "glm",
                                                predict_data = sim_data_est %>% dplyr::select(-c(Y, A)))

# Outcome regression models (Outcome probability when treatment and control)
prob_outcome1_trt1_glm <- estimate_prob_individual_model(covariates_df = sim_data_train %>% dplyr::select(-c(Y)),
                                                         label_vector = sim_data_train %>% dplyr::select(Y) %>% {.[[1]]},
                                                         params = params_glm_outcome, 
                                                         model = "glm",
                                                         predict_data = sim_data_train %>% dplyr::select(-c(Y)) %>% mutate(A = 1))
prob_outcome1_trt0_glm <- estimate_prob_individual_model(covariates_df = sim_data_train %>% dplyr::select(-c(Y)),
                                                         label_vector = sim_data_train %>% dplyr::select(Y) %>% {.[[1]]},
                                                         params = params_glm_outcome, 
                                                         model = "glm",
                                                         predict_data = sim_data_train %>% dplyr::select(-c(Y)) %>% mutate(A = 0))


estimate_AIPW_glm <- 
  estimate_AIPW(
    trt = sim_data_est$A,
    outcome = sim_data_est$Y,
    prob_trt1 = prob_trt1_glm,
    prob_outcome1_trt1 = prob_outcome1_trt1_glm,
    prob_outcome1_trt0 = prob_outcome1_trt0_glm
  )



# HOIF computation
HOIF_glm_eff1 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_glm, type = 'trt'),
  y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_glm, type = 'outcome'),
  basis = rep(basis_est,3),
  sigma = c(sigma_oracle_eff1, sigma_train_eff1, sigma_est_eff1),
  num_cores = 1
)

# Manual computation
a_resid = compute_resid(trt = sim_data_est$A, pred = prob_trt1_glm, type = 'trt')
y_resid = compute_resid(trt = sim_data_est$A, outcome = sim_data_est$Y, pred = prob_outcome1_trt1_glm, type = 'outcome')
basis_list = rep(basis_est,3)
sigma_list = c(sigma_oracle_eff1, sigma_train_eff1, sigma_est_eff1)

### Manually for one sigma - basis with one dimension - one interaction - 
basis <- basis_list[[1]]
sigma <- sigma_list[[1]]

n <- length(y_resid)
sigma_svd <- svd(sigma)
sigma_inv_sqrt <- eigenMapMatMult(eigenMapMatMult(sigma_svd$v, diag(sigma_svd$d^(-1/2))), t(sigma_svd$u))
basis_sigma_inv_sqrt <- eigenMapMatMult(basis, sigma_inv_sqrt)

y_term <- as.numeric(eigenMapMatMult(t(y_resid), basis_sigma_inv_sqrt))
a_term <- as.numeric(eigenMapMatMult(t(a_resid), basis_sigma_inv_sqrt))

stat_1 <-  ( (sum(y_term * a_term) - 
              sum(y_resid * rowSums(basis_sigma_inv_sqrt^2) * a_resid)) / (n * (n-1)) ) * (-1)

## Subtracting the diagonal terms 

### Manually for one sigma - basis with two dimension - two interaction - 
basis <- basis_list[[2]]
sigma <- sigma_list[[2]]

n <- length(y_resid)
sigma_svd <- svd(sigma)
sigma_inv_sqrt <- eigenMapMatMult(eigenMapMatMult(sigma_svd$v, diag(sigma_svd$d^(-1/2))), t(sigma_svd$u))
basis_sigma_inv_sqrt <- eigenMapMatMult(basis, sigma_inv_sqrt)

y_term <- as.numeric(eigenMapMatMult(t(y_resid), basis_sigma_inv_sqrt))
a_term <- as.numeric(eigenMapMatMult(t(a_resid), basis_sigma_inv_sqrt))

stat_2 <-  ( (sum(y_term * a_term) - 
              sum(y_resid * rowSums(basis_sigma_inv_sqrt^2) * a_resid)) / (n * (n-1)) ) * (-1)



# Compute for first sigma
HOIF_glm_eff0 <- compute_HOIF_sequence(
  a_resid = compute_resid(trt = 1-sim_data_est$X, pred = 1-prob_trt1_glm, type = 'trt'),
  y_resid = compute_resid(trt = 1-sim_data_est$X, outcome = sim_data_est$Y, pred = prob_outcome1_trt0_glm, type = 'outcome'),
  basis = rep(basis_est, 7),
  sigma = c(sigma_oracle_eff0, sigma_tr_eff0, sigma_est_eff0),
  num_cores = 1
)
