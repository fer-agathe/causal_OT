library(tidyverse)
library(ks)
library(dichromat)
library(fairadapt)
library(expm)
library(cluster)
library(transportsimplex)
library(ggtern)
library(randomForest)
library(grf)

source("functions.R")
source("utils.R")

# Data----


gen_data <- function(seed) {
  set.seed(seed)
  n_0 <- 200
  n_1 <- 400
  
  # X1 and X2 in both groups from sensitive-specific multivariate normal 
  # distributions
  M_0 <- c(-1, -1)
  S_0 <- matrix(c(1, .5, .5, 1) * 1.2^2, 2, 2)
  M_1 <- c(1.5, 1.5)
  S_1 <- matrix(c(1, -.4, -.4, 1) * 0.9^2, 2, 2)
  X_0 <- MASS::mvrnorm(n = n_0, mu = M_0, Sigma = S_0)
  X_1 <- MASS::mvrnorm(n = n_1, mu = M_1, Sigma = S_1)
  
  # Counterfactuals
  X_0_cf <- MASS::mvrnorm(n = n_0, mu = M_1, Sigma = S_1)
  X_1_cf <- MASS::mvrnorm(n = n_1, mu = M_0, Sigma = S_0)
  
  # X3: categorical, depends on S, X1, X3
  scores <- function(x1, x2, s) {
    p_A <- 0.5 + 0.3 * x1 - 0.4 * x2 + 0.2 * s
    p_B <- -0.3 + 0.5 * x2 - 0.2*x1 - 0.1 * s
    p_C <- 0
    exps <- exp(cbind(p_A, p_B, p_C))
    prob <- exps / rowSums(exps)
    prob
  }
  
  prob_X3_0 <- scores(x1 = X_0[, 1], x2 = X_0[, 2], s = 0)
  prob_X3_1 <- scores(x1 = X_1[, 1], x2 = X_1[, 2], s = 1)
  X3_0 <- apply(prob_X3_0, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  X3_1 <- apply(prob_X3_1, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  
  # Counterfactuals
  prob_X3_0_cf <- scores(x1 = X_0_cf[, 1], x2 = X_0_cf[, 2], s = 1)
  prob_X3_1_cf <- scores(x1 = X_1_cf[, 1], x2 = X_1_cf[, 2], s = 0)
  X3_0_cf <- apply(prob_X3_0_cf, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  X3_1_cf <- apply(prob_X3_1_cf, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  
  
  # Predictor for Y:
  eta_0 <- -0.2 + 0.6 * X_0[, 1] - 0.6 * X_0[, 2] + 
    ifelse(X3_0 == "B", 0.2, ifelse(X3_0 == "C", -0.3, 0))
  eta_1 <- 0.1 - 0.2 * X_1[, 1] + 0.8 * X_1[, 2] + 
    ifelse(X3_1 == "B", -0.2, ifelse(X3_1 == "C", -0.1, 0))
  
  p_0 <- exp(eta_0) / (1 + exp(eta_0))
  p_1 <- exp(eta_1) / (1 + exp(eta_1))
  
  # Predictor for Y, counterfactuals
  eta_0_cf <- 0.1 - 0.2 * X_0_cf[, 1] + 0.8 * X_0_cf[, 2] + 
    ifelse(X3_0_cf == "B", -0.2, ifelse(X3_0_cf == "C", -0.1, 0))
  
  eta_1_cf <- -0.2 + 0.6 * X_1_cf[, 1] - 0.6 * X_1_cf[, 2] + 
    ifelse(X3_1_cf == "B", 0.2, ifelse(X3_1_cf == "C", -0.3, 0))
  
  p_0_cf <- exp(eta_0_cf) / (1 + exp(eta_0_cf))
  p_1_cf <- exp(eta_1_cf) / (1 + exp(eta_1_cf))
  
  
  Y_0 <- rbinom(n_0, size = 1, prob = p_0)
  Y_1 <- rbinom(n_1, size = 1, prob = p_1)
  
  Y_0_cf <- rbinom(n_0, size = 1, prob = p_0_cf)
  Y_1_cf <- rbinom(n_1, size = 1, prob = p_1_cf)
  
  # Dataset with individuals in group 0 only
  D_SXY0 <- tibble(
    S = 0, 
    X1 = X_0[, 1], X2 = X_0[, 2], X3 = X3_0, Y = Y_0,
    X1_cf = X_0_cf[, 1], X2_cf = X_0_cf[, 2], X3_cf = X3_0_cf, Y_cf = Y_0_cf,
    eta = eta_0, p = p_0,
    eta_cf = eta_0_cf, p_cf = p_0_cf
  ) |> 
    bind_cols(as_tibble(prob_X3_0) |> rename_with(~ str_c("X3_", .))) |> 
    bind_cols(as_tibble(prob_X3_0_cf) |> rename_with(~ str_c("X3_cf_", .)))
  # Dataset with individuals in group 1 only
  D_SXY1 <- tibble(
    S = 1, 
    X1 = X_1[, 1], X2 = X_1[, 2], X3 = X3_1, Y = Y_1,
    X1_cf = X_1_cf[, 1], X2_cf = X_1_cf[, 2], X3_cf = X3_1_cf, Y_cf = Y_1_cf,
    eta = eta_1,  p = p_1,
    eta_cf = eta_1_cf, p_cf = p_1_cf
  ) |> 
    bind_cols(as_tibble(prob_X3_1) |> rename_with(~ str_c("X3_", .))) |> 
    bind_cols(as_tibble(prob_X3_1_cf) |> rename_with(~ str_c("X3_cf_", .)))
  # # Combine final dataset
  D_SXY <- rbind(D_SXY0, D_SXY1)
  
  D_SXY
}


tb <- gen_data(2)
ggplot(
  data = tb |> mutate(S = factor(S)) |> select(S, p, p_cf) |> 
    pivot_longer(cols = c(p, p_cf), names_to = "type", values_to = "p") |> 
    mutate(
      type = factor(type, levels = c("p", "p_cf"), labels = c("Obs.", "Counterfactual")
      )
    ),
  mapping = aes(x = p)
) +
  geom_histogram(
    mapping = aes(fill = S), alpha = .5, colour = "black",
    position = "identity"
  ) +
  facet_wrap(~type) +
  scale_fill_manual(values = c("0" = colours[["0"]], "1" = colours[["1"]])) +
  theme_paper()




variables <- c("S", "X1", "X2", "X3", "Y")

adj <- matrix(
  # S  X1 X2 X3 Y
  c(0, 1, 1, 1, 1,# S
    0, 0, 1, 1, 1,# X1
    0, 0, 0, 1, 1,# X2
    0, 0, 0, 0, 1,# X3
    0, 0, 0, 0, 0  # Y
  ),
  ncol = length(variables),
  dimnames = rep(list(variables), 2),
  byrow = TRUE
)

causal_graph <- fairadapt::graphModel(adj)
# plot(causal_graph)


# data <- D_SXY
# s <-  "S"
# S_0 <-  0
# y <- "Y"

# data$S |> table()

# From S=1 to S=0
library(pbapply)
library(parallel)
ncl <- detectCores()-1
(cl <- makeCluster(ncl))

clusterEvalQ(cl, {
  library(transportsimplex)
}) |>
  invisible()

seeds <- 1:200
res_simul <- vector(mode = "list", length = length(seeds))
res_ATT <- vector(mode = "list", length = length(seeds))


for (i in 1:length(seeds)) {
  cat(paste0("Simulation ", i, "/", length(seeds), "\n"))
  seed <- seeds[i]
  tb_all <- gen_data(seed)
  tb <- tb_all |> select(S, X1, X2, X3, Y) |> 
    mutate(across(is.character, ~as.factor(.x)))
  
  S_name <- "S"
  S_0 <- 0
  Y_name <- "Y"
  
  sequential_transport <- seq_trans(
    data = tb, 
    adj = adj, 
    s = S_name, 
    S_0 = S_0, 
    y = Y_name, 
    num_neighbors = 50, 
    num_neighbors_q = NULL,
    silent = FALSE,
    cl = cl
  )
  res_simul[[i]] <- sequential_transport
  
  tb_0 <- tb |> filter(!!sym(S_name) == !!S_0)
  tb_1 <- tb |> filter(!!sym(S_name) != !!S_0)
  
  n_0 <- nrow(tb_0)
  n_1 <- nrow(tb_1)
  
  ## Measuring Treatement Effect----
  mu0_model <- randomForest(
    x = tb_0 |> select(-!!Y_name, -!!S_name),
    y = factor(pull(tb_0, !!Y_name), levels = c(0,1))
  )
  
  D_1_t <- sequential_transport$transported |> 
    as_tibble() |>
    unnest_wider(where(is.list))
  pred_1_t <- predict(mu0_model, newdata = D_1_t)
  
  
  set.seed(seed)
  n <- n_0 + n_1
  n_folds <- 5 #5-fold cross-fitting
  folds <- sample(rep(1:n_folds, length.out = n))
  # Init results
  ## outcomes
  mu0_hat <- rep(NA, n)
  mu1_hat <- rep(NA, n)
  ## propensity scores
  e_hat  <- rep(NA, n)
  
  for (k in 1:n_folds) {
    idx_valid <- which(folds == k)
    idx_train <- setdiff(1:n, idx_valid)
    tb_train <- tb |> slice(idx_train)
    tb_valid <- tb |> slice(-idx_train)
    # Outcome models
    mu0_model <- randomForest(
      x = tb_train |> filter(!!sym(S_name) == !!S_0) |> 
        select(-!!Y_name, -!!S_name),
      y = tb_train |> filter(!!sym(S_name) == !!S_0) |> 
        pull(!!Y_name) |> factor(levels = c(0, 1))
    )
    mu1_model <- randomForest(
      x = tb_train |> 
        filter(!!sym(S_name) != !!S_0) |> select(-!!Y_name, -!!S_name),
      y = tb_train |> filter(!!sym(S_name) != !!S_0) |> 
        pull(!!Y_name) |> factor(levels = c(0, 1))
    )
    
    mu0_hat[idx_valid] <- predict(
      mu0_model, newdata = tb_valid |> select(-!!Y_name, -!!S_name)
    ) |> as.character() |> as.numeric()
    mu1_hat[idx_valid] <- predict(
      mu1_model, newdata = tb_valid |> select(-!!Y_name, -!!S_name)
    ) |> as.character() |> as.numeric()
    
    # Propensity model
    ps_model <- glm(
      S ~ ., data = tb_train |> select(-!!Y_name),
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
  treated_idx <- which(S != S_0)
  
  aipw_terms <- S * (Y - mu0_hat) +
    (1 - S) * (e_hat / (1 - e_hat)) * (Y - mu0_hat)
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
  ATT_dlm_se <- sqrt(var(tau_hat[treated_idx]) / sum(S == 1))
  # average_treatment_effect(fit_cf, target.sample = "treated")
  Y_1 <- tb |> filter(!!sym(S_name) != !!S_0) |> pull(!!Y_name)
  Y_1_cf <- tb_all |> filter(!!sym(S_name) != !!S_0) |> pull(Y_cf)
  
  # True ATT
  ATT_true <- mean(Y_1 - Y_1_cf)
  # ATT with the counterfactuals
  ATT_cot <- mean(Y_1 - ifelse(pred_1_t == 0, 0, 1))
  
  res_ATT[[i]] <- tibble(seed = seed, ATT_true = ATT_true, ATT_cot = ATT_cot, ATT_aipw = ATT_aipw, ATT_dml = ATT_dml)
}


save(res_ATT, res_simul, file = "res_simul.rda")

stopCluster(cl)


# 
# prop.table(table(D_SXY0$X3))
# prop.table(table(D_SXY1$X3))
# prop.table(table(sequential_transport$transported$X3))


i <- 2
ggtern(
  data = gen_data(i) |> select(S, X3_p_A,X3_p_B, X3_p_C) |> 
    mutate(S = as.character(S)) |> 
    bind_rows(
      as_tibble(res_simul[[i]]$transported_prob$X3) |> 
        rename_with(.fn = ~str_c("X3_p_", .)) |> 
        mutate(S = "1 to 0")
    ) |> 
    mutate(S = factor(S, levels = c("0", "1", "1 to 0"))),
  mapping = aes(x = X3_p_A, y = X3_p_B, z = X3_p_C)
) +
  geom_point(
    mapping = aes(colour = S), size = .5, alpha = .2
  ) +
  labs(x = "A", y = "B", z = "C") +
  scale_colour_manual(
    values = c("0" = colours[["0"]], "1" = colours[["1"]], "1 to 0" = colours[["B"]])
  ) +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_paper() +
  theme(
    legend.title = element_text(size = .8 * font_size),
    legend.text = element_text(size = .8 * font_size),
    tern.axis.vshift = .08,
    tern.axis.arrow.sep = .16,
  ) +
  theme_hidetitles() +
  guides(
    colour = guide_legend(
      override.aes = list(
        size = 1.5,
        alpha = 1
      )
    )
  )




ATT <- list_rbind(res_ATT)

ggplot(
  data = ATT |> pivot_longer(cols = -seed) |> 
    mutate(
      name = factor(
        name,
        levels = c("ATT_true", "ATT_aipw", "ATT_dml", "ATT_cot"),
        labels = c("True", "AIPW", "DML", "COT")
      )
    ),
  mapping = aes(x = value, y = name)
) +
  geom_violin(
    mapping = aes(fill = name), 
    draw_quantiles = c(0.25, 0.5, 0.75)) +
  scale_fill_manual(
    "ATT", 
    values = c(
      "True" = "#E69F00", "AIPW" = "#56B4E9", 
      "DML" = "#D55E00", "COT" = "#009E73"
    )
  ) +
  labs(y = NULL, title = "ATT (simulated data, 50 replications)") +
  theme_paper()


