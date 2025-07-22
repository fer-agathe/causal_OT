library(tidyverse)
# library(mlr3fairness)
library(seqtransfairness)
library(transportsimplex)
library(randomForest)
library(grf)
library(cluster)
source("functions.R")

data(compas, package = "mlr3fairness")


tb <- compas |> 
  as_tibble() |> 
  select(
    race, # sensitive
    age, 
    priors_count, # The prior criminal records of defendants. 
    c_charge_degree, # F: Felony M: Misdemeanor
    is_recid # outcome
  ) |> 
  mutate(
   race = ifelse(race == "Caucasian", 0, 1), # Non white as "treated"
   is_recid = ifelse(is_recid == 0, 0, 1)
  ) |> 
  sample_frac(.2)
  
  
seed <- 1234
set.seed(seed)


summary(tb)

variables <- c("race", 
               "age", "priors_count", "c_charge_degree", 
               "is_recid")

# Row: outgoing arrow
adj <- matrix(
  # S  1  2  3  Y
  c(0, 1, 1, 1, 1,# S
    0, 0, 1, 2, 1,# 1 (age)
    0, 0, 0, 0, 1,# 2 (priors_count)
    0, 0, 0, 0, 1,# 3 (c_charge_degree)
    0, 0, 0, 0, 0 # Y
  ),
  ncol = length(variables),
  dimnames = rep(list(variables), 2),
  byrow = TRUE
)

causal_graph <- fairadapt::graphModel(adj)
plot(causal_graph)

library(pbapply)
library(parallel)
ncl <- detectCores()-1
(cl <- makeCluster(ncl))

clusterEvalQ(cl, {
  library(transportsimplex)
}) |>
  invisible()


S_name <- "race"
Y_name <- "is_recid"

sequential_transport <- seq_trans(
  data = tb, 
  adj = adj, 
  s = S_name, 
  S_0 = 1, # source: treated
  y = Y_name, 
  num_neighbors = 50, 
  num_neighbors_q = NULL,
  silent = FALSE,
  cl = cl
)

stopCluster(cl)

S_untreated <- 0 # We will thus transport from 1 to 0
tb_untreated <- tb |> filter(!!sym(S_name) == !!S_untreated)
tb_treated <- tb |> filter(!!sym(S_name) != !!S_untreated)

n_untreated <- nrow(tb_untreated)
n_treated <- nrow(tb_treated)

## Measuring Treatement Effect----
mu_untreated_model <- randomForest(
  x = tb_untreated |> select(-!!Y_name, -!!S_name),
  y = factor(pull(tb_untreated, !!Y_name), levels = c(0,1))
)

D_treated_t <- sequential_transport$transported |> 
  as_tibble() |>
  unnest_wider(where(is.list))
pred_treated_t <- predict(mu_untreated_model, newdata = D_treated_t)


set.seed(seed)
n <- n_untreated + n_treated
n_folds <- 5 #5-fold cross-fitting
folds <- sample(rep(1:n_folds, length.out = n))
# Init results
## outcomes
mu_untreated_hat <- rep(NA, n)
mu_treated_hat <- rep(NA, n)
## propensity scores
e_hat  <- rep(NA, n)

for (k in 1:n_folds) {
  idx_valid <- which(folds == k)
  idx_train <- setdiff(1:n, idx_valid)
  tb_train <- tb |> slice(idx_train)
  tb_valid <- tb |> slice(-idx_train)
  # Outcome models
  mu_untreated_model <- randomForest(
    x = tb_train |> filter(!!sym(S_name) == !!S_untreated) |> 
      select(-!!Y_name, -!!S_name),
    y = tb_train |> filter(!!sym(S_name) == !!S_untreated) |> 
      pull(!!Y_name) |> factor(levels = c(0, 1))
  )
  mu_treated_model <- randomForest(
    x = tb_train |> 
    y = tb_train |> filter(!!sym(S_name) != !!S_untreated) |> 
      pull(!!Y_name) |> factor(levels = c(0, 1))
  )
  
  mu_untreated_hat[idx_valid] <- predict(
    mu_untreated_model, newdata = tb_valid |> select(-!!Y_name, -!!S_name)
  ) |> as.character() |> as.numeric()
  mu_treated_hat[idx_valid] <- predict(
    mu_treated_model, newdata = tb_valid |> select(-!!Y_name, -!!S_name)
  ) |> as.character() |> as.numeric()
  
  # Propensity model
  ps_model <- glm(
    paste(S_name, " ~ ."), data = tb_train |> select(-!!Y_name),
    family = binomial()
  )
  # Propensity scores
  e_hat[idx_valid] <- predict(
    ps_model, newdata = tb_valid, type = "response"
  )
}

# ATT with AIPW
S <- pull(tb, !!S_name)
Y <- pull(tb, !!Y_name)
treated_idx <- which(S != S_untreated)

aipw_terms <- S * (Y - mu_untreated_hat) +
  (1 - S) * (e_hat / (1 - e_hat)) * (Y - mu_untreated_hat)
ATT_aipw <- sum(aipw_terms[treated_idx]) / sum(S == 1)

# ATT with DML

# 1-hot-encoding
X_mat <- tb |> select(-!!Y_name, -!!S_name)
X_mat <- model.matrix(~ . - 1, data = X_mat)

fit_cf <- causal_forest(
  X = X_mat, 
  Y = Y, 
  W = S
)
# CATE: individual treatment effects
tau_hat <- predict(fit_cf)$predictions
ATT_dml <- mean(tau_hat[treated_idx])
ATT_dlm_se <- sqrt(var(tau_hat[treated_idx]) / n_treated)
# average_treatment_effect(fit_cf, target.sample = "treated")

# ATT with the counterfactuals
Y_treated_obs <- tb |> filter(!!sym(S_name) != !!S_untreated) |> pull(!!Y_name)
ATT_cot <- mean(Y_treated_obs - ifelse(pred_treated_t == 0, 0, 1))

tibble(ATT_cot = ATT_cot, ATT_aipw = ATT_aipw, ATT_dml = ATT_dml)
