# Application to the COMPAS recidivism dataset

# We calculate the average total causal effect, 
# the average natural direct effect,n
# and the causal mediation effect (or indirect effect) of COMPAS dataset
# assuming a known causal graph.
# We compare three methodologies to estimate those effects:
#   
# 1. Using Optimal Transport (OT) to first derive counterfactuals at the 
#    individual level and then, averaging over all the individuals in the dataset,
# 2. Using Sequential Transport (ST) to first derive counterfactuals at the 
#    individual level and then, averaging over all the individuals in the dataset,
# 3. Using LSEM from causal mediation analysis that fits linear models to 
#    estimate the different average causal effects.

# 1 Setup----

## 1.1 Some Packages----

library(tidyverse)
# install.packages("mlr3fairness")
# remotes::install_github(
#   repo = "fer-agathe/sequential_transport", subdir = "seqtransfairness"
# )
library(seqtransfairness)
# remotes::install_github(repo = "fer-agathe/transport-simplex")
library(transportsimplex)
library(randomForest)
library(fairadapt)
library(grf)
library(cluster)


## 1.2 Graphs----

library(extrafont, quietly = TRUE)
font_family <- "serif"

path <- "./figs/"
if (!dir.exists(path)) dir.create(path)

library(tikzDevice)
# Functions to export graphs
source("../scripts/utils.R")

## 1.3 Functions----

# Functions used to build counterfactuals
source("../scripts/functions.R")

# 2 Data----

data(compas, package = "fairadapt")

# The outcome variable is binary here: 
# whether the defendant is rearrested at any time. 
# The "treatment" A will be the sensitive attribute, `race`. 
# In the data compiled in {fairadapt}, race is defined as `White` 
# (which we will consider as A=1),
# and `Non-White` (which we will consider as A=0). 
# The idea is to build counterfactuals for non-white people to ask questions 
# such as "had this non-white individual been white, what would the 
# prediction of an algorithm modeling recidivism be?".

# To train the predictive model of recidivism, we will use the following 
# covariates: the age, the prior criminal records of defendants, 
# and the charge degree (felony or misdemeanor).
tb <- compas |> 
  as_tibble() |> 
  select(
    race, # sensitive
    age, 
    priors_count, # The prior criminal records of defendants. 
    c_charge_degree, # F: Felony M: Misdemeanor
    is_recid = two_year_recid # outcome
  ) |> 
  mutate(
    race = ifelse(race == "Non-White", 0, 1), # non-white indiv. as "untreated"
  )

dim(tb)
summary(tb)

# Assumed structural relationship
variables <- c("race", 
               "age", "priors_count", "c_charge_degree", 
               "is_recid")
# Row: outgoing arrow
adj <- matrix(
  # S  1  2  3  Y
  c(0, 1, 1, 1, 1,# S
    0, 0, 1, 1, 1,# 1 (age)
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

# 3 Counterfactuals----

seed <- 1234
set.seed(seed)

A_name <- "race" # treatment name
Y_name <- "is_recid" # outcome name
A_untreated <- 0
A <- tb[[A_name]]
ind_untreated <- which(A == A_untreated)
tb_untreated <- tb[ind_untreated, ]
tb_treated <- tb[-ind_untreated, ]

## 3.1 Multivariate Optimal Transport----

tb_untreated_wo_A <- tb_untreated[ , !(names(tb_untreated) %in% A_name)]
tb_treated_wo_A <- tb_treated[ , !(names(tb_treated) %in% A_name)]
n0 <- nrow(tb_untreated_wo_A)
n1 <- nrow(tb_treated_wo_A)
y0 <- tb_untreated_wo_A[[Y_name]]
y1 <- tb_treated_wo_A[[Y_name]]
X0 <- tb_untreated_wo_A[ , !(names(tb_untreated_wo_A) %in% Y_name)]
X1 <- tb_treated_wo_A[ , !(names(tb_treated_wo_A) %in% Y_name)]

# One hot encoding of categorical variables
num_cols <- names(X0)[sapply(X0, is.numeric)]
cat_cols <- names(X0)[sapply(X0, function(col) is.factor(col) || is.character(col))]
X0_num <- X0[ , num_cols]
X1_num <- X1[ , num_cols]
X0_cat <- X0[ , cat_cols]
X1_cat <- X1[ , cat_cols]

cat_counts <- sapply(X0[ , cat_cols], function(col) length(unique(col)))

library(caret)
X0_cat_encoded <- list()
X1_cat_encoded <- list()
for (col in cat_cols) {
  # One-hot encoding with dummyVars
  formula <- as.formula(paste("~", col))
  dummies <- dummyVars(formula, data = X0_cat)
  
  # Dummy variable
  dummy_0 <- predict(dummies, newdata = X0_cat) %>% as.data.frame()
  dummy_1 <- predict(dummies, newdata = X1_cat) %>% as.data.frame()
  
  # Scaling
  dummy_0_scaled <- scale(dummy_0)
  dummy_1_scaled <- scale(dummy_1)
  
  dummy_0_df <- as.data.frame(dummy_0_scaled)
  dummy_1_df <- as.data.frame(dummy_1_scaled)
  
  # Aling categories in both treated/untreated groups
  all_cols <- union(colnames(dummy_0_df), colnames(dummy_1_df))
  dummy_0_df <- dummy_0_df %>% mutate(across(.fns = identity)) %>% select(all_of(all_cols)) %>% replace(is.na(.), 0)
  dummy_1_df <- dummy_1_df %>% mutate(across(.fns = identity)) %>% select(all_of(all_cols)) %>% replace(is.na(.), 0)
  
  # Sauvegarde dans les listes
  X0_cat_encoded[[col]] <- dummy_0_df
  X1_cat_encoded[[col]] <- dummy_1_df
}

# Euclidean distance for numerical variables
num_dist <- proxy::dist(x = X0_num, y = X1_num, method = "Euclidean")
num_dist <- as.matrix(num_dist)

# Hamming distance for categorical variables
cat_dists <- list()
for (col in cat_cols) {
  mat_0 <- as.matrix(X0_cat_encoded[[col]])
  mat_1 <- as.matrix(X1_cat_encoded[[col]])
  dist_mat <- proxy::dist(x = mat_0, y = mat_1, method = "Euclidean")
  cat_dists[[col]] <- as.matrix(dist_mat)
}

# Combine the two distance matrices
combined_cost <- num_dist
for (i in seq_along(cat_dists)) {
  combined_cost <- combined_cost + cat_dists[[i]]
}

# Transport map:
# Uniform weights (equal mass)
w0 <- rep(1 / n0, n0)
w1 <- rep(1 / n1, n1)
# Compute transport plan
transport_res <- transport::transport(
  a = w0,
  b = w1,
  costm = combined_cost,
  method = "shortsimplex"
)

# Solving the transport problem
transport_plan <- matrix(0, nrow = n0, ncol = n1)
for(i in seq_len(nrow(transport_res))) {
  transport_plan[transport_res$from[i], transport_res$to[i]] <- transport_res$mass[i]
}

# Using the transport plan for numerical variables:
num_transported <- n0 * (transport_plan %*% as.matrix(X1_num))

# and for categorical variables (+label reconstruction)
cat_transported <- list()
for (col in cat_cols) {
  cat_probs <- transport_plan %*% as.matrix(X1_cat_encoded[[col]])
  cat_encoded_columns <- colnames(X1_cat_encoded[[col]])
  # For each obs., we take the index with the maximum value (approx. proba)
  max_indices <- apply(cat_probs, 1, which.max)
  prefix_pattern <- paste0("^", col, "\\.")
  cat_transported[[col]] <- sapply(max_indices, function(x) sub(prefix_pattern, "", cat_encoded_columns[x]))
}

# Put the result in a tibble
tb_ot_transported <- as_tibble(num_transported)
for (col in cat_cols) {
  tb_ot_transported[[col]] <- cat_transported[[col]]
}

tb_ot_transported <- tb_ot_transported |>
  mutate(c_charge_degree = as.factor(c_charge_degree))
tb_ot_transported <- as.list(tb_ot_transported)

## 3.2 Penalized Optimal Transport----

# Transport map using Sinkhorn penalty
sinkhorn_transport_res <- T4transport::sinkhornD(
  combined_cost, wx = w0, wy = w1, lambda = 0.1
)
# The estimated transport plan
sinkhorn_transport_plan <- sinkhorn_transport_res$plan

# Transport numeric variables
num_sinkhorn_transported <- n0 * (sinkhorn_transport_plan %*% as.matrix(X1_num))

# Then categorical variables
cat_sinkhorn_transported <- list()
for (col in cat_cols) {
  cat_probs <- sinkhorn_transport_plan %*% as.matrix(X1_cat_encoded[[col]])
  cat_encoded_columns <- colnames(X1_cat_encoded[[col]])
  # For each obs., we take the index with the maximum value (approx. proba)
  max_indices <- apply(cat_probs, 1, which.max)
  prefix_pattern <- paste0("^", col, "\\.")
  cat_sinkhorn_transported[[col]] <- sapply(max_indices, function(x) sub(prefix_pattern, "", cat_encoded_columns[x]))
}

# Store results in a tibble
tb_sinkhorn_transported <- as_tibble(num_sinkhorn_transported)
for (col in cat_cols) {
  tb_sinkhorn_transported[[col]] <- cat_sinkhorn_transported[[col]]
}

tb_sinkhorn_transported <- tb_sinkhorn_transported |>
  mutate(c_charge_degree = as.factor(c_charge_degree))
tb_sinkhorn_transported <- as.list(tb_sinkhorn_transported)


## 3.3 Sequential Transport----

library(pbapply)
library(parallel)
ncl <- detectCores()-1
(cl <- makeCluster(ncl))

clusterEvalQ(cl, {
  library(transportsimplex)
}) |>
  invisible()


A_name <- "race" # treatment name
Y_name <- "is_recid" # outcome name

sequential_transport <- seq_trans(
  data = tb, 
  adj = adj, 
  s = A_name, 
  S_0 = 0, # source: untreated
  y = Y_name, 
  num_neighbors = 50, 
  num_neighbors_q = NULL,
  silent = FALSE,
  cl = cl
)

stopCluster(cl)

## 3.4 Fairadapt----

# RF
fpt_model_rf <- fairadapt(
  priors_count ~ ., 
  train.data = tb,
  prot.attr = A_name, adj.mat = adj,
  quant.method = rangerQuants
)

adapt_df_rf <- adaptedData(fpt_model_rf)
adapt_df_rf_untreated <- adapt_df_rf[ind_untreated, ]

# Linear model
fpt_model_lin <- fairadapt(
  priors_count ~ ., 
  train.data = tb,
  prot.attr = A_name, adj.mat = adj,
  quant.method = linearQuants
)

adapt_df_lin <- adaptedData(fpt_model_lin)
adapt_df_lin_untreated <- adapt_df_lin[ind_untreated, ]

# 4 Measuring the Causal Effect----

## 4.1 With Causal Mediation Analysis----
# library(mediation) # we do not load it
#  otherwise it masks a lot of useful functions

# We encode the categorical variable as for optimal transport
tb_med <- tb |> 
  mutate(
    c_charge_degree = ifelse(c_charge_degree == "F", 0, 1)
  )

med_mod_age <- mediation::multimed(
  outcome = "is_recid", 
  med.main = "age", 
  med.alt = c("priors_count", "c_charge_degree"), 
  treat = "race", 
  data = tb_med
)
# Indirect effect for age: race -> age -> ir
delta_0_med_age <- mean((med_mod_age$d0.lb + med_mod_age$d0.ub) / 2)
# Direct + Other indirect effects: race -> ir, race -> pc -> ir, race -> ccd -> ir, 
# race -> age -> pc -> ir, race -> age -> ccd -> ir
zeta_1_med_age <- mean((med_mod_age$z1.lb + med_mod_age$z1.ub) / 2)
# Total effect
tot_effect_med_age <- delta_0_med_age + zeta_1_med_age

med_mod_pc <- mediation::multimed(
  outcome = "is_recid", 
  med.main = "priors_count", 
  med.alt = c("age", "c_charge_degree"), 
  treat = "race", 
  data = tb_med
)
# Indirect effect for pc: race -> pc -> ir, race -> age -> pc -> ir
delta_0_med_pc <- mean((med_mod_pc$d0.lb + med_mod_pc$d0.ub) / 2)
# Direct + Other indirect effects: race -> ir, race -> age -> ir, race -> ccd -> ir, 
# race -> age -> ccd -> ir
zeta_1_med_pc <- mean((med_mod_pc$z1.lb + med_mod_pc$z1.ub) / 2)
# Total effect
tot_effect_med_pc <- delta_0_med_pc + zeta_1_med_pc

med_mod_ccd <- mediation::multimed(
  outcome = "is_recid", 
  med.main = "c_charge_degree", 
  med.alt = c("age", "priors_count"), 
  treat = "race", 
  data = tb_med
)
# Indirect effect for ccd: race -> ccd -> ir, race -> age -> ccd -> ir
delta_0_med_ccd <- mean((med_mod_ccd$d0.lb + med_mod_ccd$d0.ub) / 2)
# Direct + Other indirect effects: race -> ir, race -> age -> ir, race -> pc -> ir, 
# race -> age -> pc -> ir
zeta_1_med_ccd <- mean((med_mod_ccd$z1.lb + med_mod_ccd$z1.ub) / 2)
# Total effect
tot_effect_med_ccd <- delta_0_med_ccd + zeta_1_med_ccd

# The estimated values:
# Total effect
tot_effect_med <- tot_effect_med_age
# Indirect effects
delta_0_med <- delta_0_med_age + delta_0_med_pc + delta_0_med_ccd
# Direct effect 
zeta_1_med <- tot_effect_med - delta_0_med
cbind(delta_0 = delta_0_med, zeta_1 = zeta_1_med, tot_effect = tot_effect_med)

## 4.2 With Optimal Transport----

library(randomForest)
# Random forest to estimate the outcome model
causal_effects_ot <- causal_effects_cf(
  data_untreated = tb_untreated, 
  data_treated = tb_treated,
  data_cf_untreated = as_tibble(tb_ot_transported), 
  data_cf_treated = NULL, 
  Y_name = Y_name, 
  A_name = A_name, 
  A_untreated = A_untreated # 0
)

cbind(
  delta_0 = causal_effects_ot$delta_0,
  zeta_1 = causal_effects_ot$zeta_1, 
  tot_effect = causal_effects_ot$tot_effect
)

## 4.3 With Penalized Optimal Transport----
causal_effects_sink_ot <- causal_effects_cf(
  data_untreated = tb_untreated, 
  data_treated = tb_treated,
  data_cf_untreated = as_tibble(tb_sinkhorn_transported), 
  Y_name = Y_name, 
  A_name = A_name, 
  A_untreated = A_untreated # 0
)

cbind(
  delta_0 = causal_effects_sink_ot$delta_0,
  zeta_1 = causal_effects_sink_ot$zeta_1, 
  tot_effect = causal_effects_sink_ot$tot_effect
)

## 4.4 With Sequential Transport----
causal_effects_st <- causal_effects_cf(
  data_untreated = tb_untreated, 
  data_treated = tb_treated,
  data_cf_untreated = as_tibble(sequential_transport$transported), 
  Y_name = Y_name, 
  A_name = A_name, 
  A_untreated = A_untreated
)

cbind(
  delta_0 = causal_effects_st$delta_0,
  zeta_1 = causal_effects_st$zeta_1, 
  tot_effect = causal_effects_st$tot_effect
)

## 4.5 With fairadapt----
causal_effects_fairadapt_rf <- causal_effects_cf(
  data_untreated = tb_untreated, 
  data_treated = tb_treated,
  data_cf_untreated = as_tibble(adapt_df_rf_untreated |> select(-c(!!A_name, !!Y_name))), 
  Y_name = Y_name, 
  A_name = A_name, 
  A_untreated = A_untreated
)

causal_effects_fairadapt_lin <- causal_effects_cf(
  data_untreated = tb_untreated, 
  data_treated = tb_treated,
  data_cf_untreated = as_tibble(adapt_df_lin_untreated |> select(-c(!!A_name, !!Y_name))), 
  Y_name = Y_name, 
  A_name = A_name, 
  A_untreated = A_untreated
)

cbind(
  delta_0_rf = causal_effects_fairadapt_rf$delta_0,
  zeta_1_rf = causal_effects_fairadapt_rf$zeta_1, 
  tot_effect_rf = causal_effects_fairadapt_rf$tot_effect,
  delta_0_lin = causal_effects_fairadapt_lin$delta_0,
  zeta_1_lin = causal_effects_fairadapt_lin$zeta_1, 
  tot_effect_lin = causal_effects_fairadapt_lin$tot_effect
)

## 4.6 With Normalizing Flows----

# First export the dataset
write_csv(tb, file = "../output/compas.csv")

# Then create a python environment and install required python libraries
# This is done in python (see script `NF_compas.py`)
system("python3 -m venv venv")
system("venv/bin/python -m pip install --upgrade pip")
system(paste0(
  "venv/bin/python -m pip install ",
  "--extra-index-url https://test.pypi.org/simple/ medflow"
))
system("venv/bin/python -m pip install numpy")

# Execute the script (this takes about 30 min on a standard computer, in 2026)
# If this code is executed in a terminal instead of called from R with system()
# a progress is shown.
system("venv/bin/python 03_NF_compas.py")

# We load the obtained results:
nf_estim <- read.csv("../output/1_path_100k_v2_results.csv")

nf_tot_effect <- 
  nf_estim |> filter(Potential.Outcome == "E[is_recid(race=1)]") |> pull("Value") -
  nf_estim |> filter(Potential.Outcome == "E[is_recid(race=0)]") |> pull("Value")

# Indirect
nf_delta_0 <- 
  nf_estim |> filter(
    Potential.Outcome == "E[is_recid(race=0, age(race=1), priors_count(race=1), c_charge_degree(race=1))]"
  ) |> pull("Value") -
  nf_estim |> filter(
    Potential.Outcome == "E[is_recid(race=0)]"
  ) |> pull("Value")

# Direct
nf_zeta_1 <- nf_tot_effect - nf_delta_0

# Sanity check for direct effect
nf_zeta_1_check <- nf_estim |> filter(
  Potential.Outcome == "E[is_recid(race=1)]"
) |> pull("Value") -
  nf_estim |> filter(
    Potential.Outcome == "E[is_recid(race=0, age(race=1), priors_count(race=1), c_charge_degree(race=1))]"
  ) |> pull("Value")

c(nf_zeta_1, nf_zeta_1_check)


# 5 Visualization of the Results----

## 5.1 Histograms for individual effects (Figure 4)----

# Note that the individual effects are not provided for NF.

plot_hist_effects <- function(x, 
                              var_name, 
                              tikz = FALSE,
                              fill = "red",
                              printed_method = "",
                              x_lim = NULL,
                              print_main = TRUE,
                              print_x_axis = TRUE) {
  if (print_main == TRUE) {
    name_effect <- case_when(
      str_detect(var_name, "^delta_0") ~ "$\\delta_i$",#"$\\delta_i(0)$",
      str_detect(var_name, "^zeta_1") ~ "$\\zeta_i$",#"$\\zeta_i(1)$",
      str_detect(var_name, "^tot_effect") ~ "$\\tau_i$",#"$\\tau_i(1)$",
      TRUE ~ "other"
    )
    if (tikz == FALSE) name_effect <- latex2exp::TeX(name_effect)
  } else {
    name_effect <- ""
  }
  
  
  if (var_name == "tot_effect") {
    data_plot <- x[["delta_0_i"]] + x[["zeta_1_i"]]
  } else {
    data_plot <- x[[var_name]]
  }
  
  if (is.null(x_lim)) {
    hist(
      data_plot, 
      main = "", xlab = "", ylab = "", family = font_family,
      col = fill, axes = FALSE
    )
  } else {
    hist(
      data_plot, 
      main = "", xlab = "", ylab = "", family = font_family,
      col = fill, xlim = x_lim, axes = FALSE
    )
  }
  
  if (print_x_axis) axis(1, family = font_family)
  axis(2, family = font_family)
  
  title(
    main = name_effect, cex.main = 1, family = font_family
  )
  
  if (printed_method != "") {
    title(
      ylab = printed_method, line = 2, 
      cex.lab = 1, family = font_family
    )
  }
  abline(v = mean(data_plot), col = "darkred", lty = 2, lwd = 2)
}


colour_methods <- c(
  "fairadapt_rf" = "#9966FF",
  # "OT" = "#CC79A7",
  "OT-M" = "#009E73",
  "skh" = "darkgray",
  "seq_1" = "#0072B2"
  #"seq_2" = "#D55E00",
)


export_tikz <- FALSE

scale <- 1.42
file_name <- "compas-dist-indiv-effects"
width_tikz <- 2.25*scale
height_tikz <- 1.4*scale
if (export_tikz == TRUE)
  tikz(paste0("figs/", file_name, ".tex"), width = width_tikz, height = height_tikz)

layout(
  matrix(1:12, byrow = TRUE, ncol = 3),
  widths = c(1, rep(.9, 2)), heights = c(1, rep(.72, 2))
)

for (i in 1:4) {
  x <- case_when(
    i == 1 ~ causal_effects_fairadapt_rf,
    i == 2 ~ causal_effects_ot,
    i == 3 ~ causal_effects_sink_ot,
    i == 4 ~ causal_effects_st
    
  )
  method <- case_when(
    i == 1 ~ "FPT-RF",
    i == 2 ~ "OT-M",
    i == 3 ~ "SKH",
    i == 4 ~ "ST",
  )
  
  for (var_name in c("delta_0_i", "zeta_1_i", "tot_effect")) {
    mar_bottom <- ifelse(i == 4, 2.1, .6)
    mar_left <- ifelse(var_name == "delta_0_i", 3.1, 2.1)
    mar_top <- ifelse(i == 1, 2.1, .1)
    mar_right <- .4
    printed_method <- ifelse(
      var_name == "delta_0_i", 
      ifelse(method == "OT-M", "OT", method),
      ""
    )
    
    par(mar = c(mar_bottom, mar_left, mar_top, mar_right))
    x_lim_list <- list(
      "delta_0_i" = c(-.5, .5),
      "zeta_1_i" = c(-.4, .4),
      "tot_effect" = c(-.6, .6)
    )
    
    plot_hist_effects(
      x = x, var_name = var_name, tikz = export_tikz, 
      fill = colour_methods[i],
      printed_method = printed_method, 
      x_lim = x_lim_list[[var_name]],
      print_main = i == 1,
      print_x_axis = i == 4
    )
  }
}


if (export_tikz == TRUE) {
  dev.off()
  plot_to_pdf(
    filename = file_name, 
    path = "./figs/", keep_tex = FALSE, crop = T
  )
}

## 5.2 Summary Table (Table 1)----


tibble::tribble(
  ~method, ~delta_0, ~zeta_1, ~tot_effect,
  "CM", delta_0_med, zeta_1_med, tot_effect_med,
  "FPT-RF", causal_effects_fairadapt_rf$delta_0, causal_effects_fairadapt_rf$zeta_1, causal_effects_fairadapt_rf$tot_effect,
  "FPT-LIN", causal_effects_fairadapt_lin$delta_0, causal_effects_fairadapt_lin$zeta_1, causal_effects_fairadapt_lin$tot_effect,
  "NF", nf_delta_0, nf_zeta_1, nf_tot_effect,
  "OT", causal_effects_ot$delta_0, causal_effects_ot$zeta_1, causal_effects_ot$tot_effect,
  "SKH", causal_effects_sink_ot$delta_0, causal_effects_sink_ot$zeta_1, causal_effects_sink_ot$tot_effect,
  "ST", causal_effects_st$delta_0, causal_effects_st$zeta_1, causal_effects_st$tot_effect
) |> 
  knitr::kable(digits = 3)


## 5.3 Decomposition of the Indirect Effect----

# Total indirect effect with the Sequential Transport approach:
causal_effects_st$delta_0

# We need to retrieve the fitted RF model of the total indirect effect.
data_untreated <- tb_untreated
# All variables transported
data_cf_untreated <- as_tibble(sequential_transport$transported)
Y_name <- Y_name
A_name <- A_name
# Outcome model for untreated
mu_untreated_model <- randomForest(
  x = data_untreated |> dplyr::select(-!!Y_name, -!!A_name),
  y = pull(data_untreated, !!Y_name)
)

# Now we can modify the causal_effects_cf() function to take the fitted RF 
# model as an argument to be able to decompose the indirect effect and we 
# return only the indirect effect at the individual level, `delta_0_i` and the 
# average `delta_0`.

modified_causal_effects_cf <- function(data_untreated,
                                       data_cf_untreated,
                                       mu_untreated_model) {
  
  n_untreated <- nrow(data_untreated)
  
  # Natural Indirect Effect, using predictions
  delta_0_i <- predict(mu_untreated_model, newdata = data_cf_untreated) -
    predict(mu_untreated_model, newdata = data_untreated)
  delta_0 <- mean(delta_0_i)
  
  list(
    delta_0_i = delta_0_i,
    delta_0 = delta_0
  )
}

# Total indirect effect
modified_causal_effects_cf(data_untreated, data_cf_untreated, mu_untreated_model)$delta_0

### Influence of Age----
cf_treated_first <- data_untreated |> 
  mutate(age = as_tibble(sequential_transport$transported)$age)

indirect_age <- modified_causal_effects_cf(
  data_untreated = data_untreated, 
  data_cf_untreated = cf_treated_first, 
  mu_untreated_model = mu_untreated_model
)

cbind(
  delta_0_age = indirect_age$delta_0
)


### Influence of Prior Counts----
cf_treated_second <- cf_treated_first |> 
  mutate(priors_count = as_tibble(sequential_transport$transported)$priors_count)

indirect_priors_count <- modified_causal_effects_cf(
  data_untreated = cf_treated_first, 
  data_cf_untreated = cf_treated_second, 
  mu_untreated_model = mu_untreated_model
)

cbind(
  delta_0_priors_count = indirect_priors_count$delta_0
)

### Influence of Charge Degree----

cf_treated_third <- cf_treated_second |> 
  mutate(c_charge_degree = as_tibble(sequential_transport$transported)$c_charge_degree)

indirect_c_charge_degree <- modified_causal_effects_cf(
  data_untreated = cf_treated_second, 
  data_cf_untreated = cf_treated_third, 
  mu_untreated_model = mu_untreated_model
)

cbind(
  delta_0_c_charge_degree = indirect_c_charge_degree$delta_0
)

