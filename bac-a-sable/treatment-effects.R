library(dplyr)
set.seed(1234)

# Data----

n = 1000
X1 = rnorm(n)
X2star = 1 + X1 + rnorm(n)
X2 = cut(X2star, c(-Inf, -1, +1, +Inf), c("A", "B", "C"))
X3 = X1 + (X2=="A")+rnorm(n)
eta = X1+(X2=="B") - 0.4
probaT = exp(eta)/(1+exp(eta))
T = rbinom(n, size=1, prob = probaT)
Y = 1 + 2*X1 + rnorm(n) + T*(X3+(X1=="C"))
base = tibble(
  X1 = X1,
  X2 = as.factor(X2),
  X3 = X3,
  T  = T,
  Y = Y
)

# Train/Test splits
train_idx <- sample(seq_len(n), size = floor(0.8 * n))
base_train <- base[train_idx, ]
base_test  <- base[-train_idx, ]

# One-hot encode
X_train <- model.matrix(~ X1 + X2 + X3 - 1, data = base_train)
X_test  <- model.matrix(~ X1 + X2 + X3 - 1, data = base_test)

# ATE----

## 1. Using OLS----
lm_ate <- lm(Y ~ T, data = base_train)
ate_lm <- coef(lm_ate)[["T"]]

## 2. Using a causal forest----

# The Causal Forest use honnest estimation by default...

library(grf)
fit_cf <- causal_forest(
  X = X_train, 
  Y = base_train$Y, 
  W = base_train$T
)
ate_cf <- average_treatment_effect(fit_cf, target.sample = "all")


# CATE----

## 1. X-Learner (with RF)----

# Estimation of mu_0 and mu_1
library(randomForest)
ind_train_0 <- base_train$T == 0
ind_train_1 <- base_train$T == 1
ind_test_0 <- base_test$T == 0
ind_test_1 <- base_test$T == 1

mu1_model <- randomForest(
  x = X_train[ind_train_1, ],
  y = base_train$Y[ind_train_1]
)
mu0_model <- randomForest(
  x = X_train[ind_train_0, ],
  y = base_train$Y[ind_train_0]
)

# Treatment effects:
# Treated: T_i = Y_i - \mu_0(X_i)
# mu0_hat_test <- predict(mu0_model, newdata = X_test) # potential outcome on test
# pseudo treatment effects
D1 <- base_train$Y[ind_train_1] - 
  predict(mu0_model, newdata = X_train[ind_train_1, ])

# Control: T_i = \mu_1(X_i) - Y_i
# mu1_hat_test <- predict(mu1_model, newdata = X_test) # potential outcome on test
D0 <- predict(mu1_model, newdata = X_train[ind_train_0, ]) -
  base_train$Y[ind_train_0]

# Train tau models:
tau1_model <- randomForest(x = X_train[ind_train_1, ], y = D1)
tau0_model <- randomForest(x = X_train[ind_train_0, ], y = D0)

# Propensity scores (on training set):
ps_model <- glm(T ~ X1 + X2 + X3, data = base_train, family = binomial)
g_hat_test <- predict(ps_model, newdata = base_test, type = "response")

# Ppseudo-effects on test set:
tau1_hat <- predict(tau1_model, newdata = X_test)
tau0_hat <- predict(tau0_model, newdata = X_test)

cate_xlearner <- g_hat_test * tau0_hat + (1 - g_hat_test) * tau1_hat

## 2. Causal Forest----
pred_cf <- predict(fit_cf)
cate_cf <- pred_cf$predictions


# ATT----

## 1. X_Learner (with RF)----
att_xlearner <- mean(cate_xlearner[ind_test_1])


## 2. Using Prop. Score Weighting (IPW)----

ps_model <- glm(T ~ X1 + X2 + X3, data = base_train, family = binomial)

att_ipw <- base_test |> 
  mutate(
    # Propensity score
    ps = predict(ps_model, newdata = base_test, type = "response"),
    # Weights: ps/(1-ps) for control
    weight = if_else(condition = T == 1, true = 1, false = ps / (1 - ps))
  ) |> 
  summarise(
    treated_mean = mean(Y[T == 1]),
    control_weighted_mean = weighted.mean(Y[T == 0], w = weight[T == 0]),
    ATT = treated_mean - control_weighted_mean
  ) |> 
  pull("ATT")

## 3. Causal Forest----
att_cf <- average_treatment_effect(fit_cf, target.sample = "treated")


