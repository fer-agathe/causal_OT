library(tidyverse)
# remotes::install_github(
#   repo = "fer-agathe/sequential_transport", subdir = "seqtransfairness"
# )
library(seqtransfairness)
# remotes::install_github(repo = "fer-agathe/transport-simplex")
library(transportsimplex)
library(randomForest)
library(grf)
library(cluster)
library(dplyr)

# Also required:
# install.packages(mlr3fairness)



library(extrafont, quietly = TRUE)
col_group <- c("#00A08A","#F2AD00", "#1b95e0")
colour_methods <- c(
    "OT" = "#CC79A7", 
    "OT-M" = "#009E73",
    "skh" = "darkgray",
    "seq_1" = "#0072B2", "seq_2" = "#D55E00"
  )
colGpe1 <- col_group[2]
colGpe0 <- col_group[1]
colGpet <- col_group[3]
loadfonts(device = "pdf", quiet = TRUE)
font_size <- 20
font_family <- "serif"

path <- "./figs/"
if (!dir.exists(path)) dir.create(path)

source("../scripts/utils.R")

source("../scripts/functions.R")
source("../scripts/utils.R")

gen_data <- function(seed) {
  set.seed(seed)
  n_0 <- 400
  n_1 <- 200
  
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
  scores <- function(x1, x2, a) {
    p_A <- 0.5 + 0.3 * x1 - 0.4 * x2 + 0.2 * a
    p_B <- -0.3 + 0.5 * x2 - 0.2*x1 - 0.1 * a
    p_C <- 0
    exps <- exp(cbind(p_A, p_B, p_C))
    prob <- exps / rowSums(exps)
    prob
  }
  
  prob_X3_0 <- scores(x1 = X_0[, 1], x2 = X_0[, 2], a = 0)
  prob_X3_1 <- scores(x1 = X_1[, 1], x2 = X_1[, 2], a = 1)
  X3_0 <- apply(prob_X3_0, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  X3_1 <- apply(prob_X3_1, 1, function(p) sample(c("A", "B", "C"), 1, prob = p))
  
  # Counterfactuals
  prob_X3_0_cf <- scores(x1 = X_0_cf[, 1], x2 = X_0_cf[, 2], a = 1)
  prob_X3_1_cf <- scores(x1 = X_1_cf[, 1], x2 = X_1_cf[, 2], a = 0)
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
  data_0 <- tibble(
    A = 0, 
    X1 = X_0[, 1], X2 = X_0[, 2], X3 = X3_0, Y = Y_0,
    X1_cf = X_0_cf[, 1], X2_cf = X_0_cf[, 2], X3_cf = X3_0_cf, Y_cf = Y_0_cf,
    eta = eta_0, p = p_0,
    eta_cf = eta_0_cf, p_cf = p_0_cf
  ) |> 
    bind_cols(as_tibble(prob_X3_0) |> rename_with(~ str_c("X3_", .))) |> 
    bind_cols(as_tibble(prob_X3_0_cf) |> rename_with(~ str_c("X3_cf_", .)))
  # Dataset with individuals in group 1 only
  data_1 <- tibble(
    A = 1, 
    X1 = X_1[, 1], X2 = X_1[, 2], X3 = X3_1, Y = Y_1,
    X1_cf = X_1_cf[, 1], X2_cf = X_1_cf[, 2], X3_cf = X3_1_cf, Y_cf = Y_1_cf,
    eta = eta_1,  p = p_1,
    eta_cf = eta_1_cf, p_cf = p_1_cf
  ) |> 
    bind_cols(as_tibble(prob_X3_1) |> rename_with(~ str_c("X3_", .))) |> 
    bind_cols(as_tibble(prob_X3_1_cf) |> rename_with(~ str_c("X3_cf_", .)))
  # # Combine final dataset
  data_all <- rbind(data_0, data_1)
  
  data_all
}

tb <- gen_data(2) |> 
  mutate(across(where(is.character), ~as.factor(.x)))

ggplot(
  data = tb |> mutate(A = factor(A)) |> dplyr::select(A, p, p_cf) |> 
    pivot_longer(cols = c(p, p_cf), names_to = "type", values_to = "p") |> 
    mutate(
      type = factor(type, levels = c("p", "p_cf"), labels = c("Obs.", "Counterfactual")
      )
    ),
  mapping = aes(x = p)
) +
  geom_histogram(
    mapping = aes(fill = A), alpha = .5, colour = "black",
    position = "identity", bins = 30
  ) +
  facet_wrap(~type) +
  scale_fill_manual(values = c("0" = colours[["0"]], "1" = colours[["1"]])) +
  theme_paper()

variables <- c("A", "X1", "X2", "X3", "Y")

adj <- matrix(
  # A  X1 X2 X3 Y
  c(0, 1, 1, 1, 1,# A
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
plot(causal_graph)

seed <- 1234
set.seed(seed)

A_name <- "A" # treatment name
Y_name <- "Y" # outcome name
A_untreated <- 0
A <- tb[[A_name]]
ind_untreated <- which(A == A_untreated)
tb_estim <- tb |> dplyr::select(A, X1, X2, X3, Y)
tb_untreated <- tb_estim[ind_untreated, ]
tb_treated <- tb_estim[-ind_untreated, ]

tb_untreated_wo_A <- tb_untreated[ , !(names(tb_untreated) %in% A_name)]
tb_treated_wo_A <- tb_treated[ , !(names(tb_treated) %in% A_name)]
n0 <- nrow(tb_untreated_wo_A)
n1 <- nrow(tb_treated_wo_A)
y0 <- tb_untreated_wo_A[[Y_name]]
y1 <- tb_treated_wo_A[[Y_name]]
X0 <- tb_untreated_wo_A[ , !(names(tb_untreated_wo_A) %in% Y_name)]
X1 <- tb_treated_wo_A[ , !(names(tb_treated_wo_A) %in% Y_name)]

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
  dummies <- caret::dummyVars(formula, data = X0_cat)
  
  # Dummy variable
  dummy_0 <- predict(dummies, newdata = X0_cat) |> as.data.frame()
  dummy_1 <- predict(dummies, newdata = X1_cat) |> as.data.frame()
  
  # Scaling
  dummy_0_scaled <- scale(dummy_0)
  dummy_1_scaled <- scale(dummy_1)

  dummy_0_df <- as.data.frame(dummy_0_scaled)
  dummy_1_df <- as.data.frame(dummy_1_scaled)
  
  # Aling categories in both treated/untreated groups
  all_cols <- union(colnames(dummy_0_df), colnames(dummy_1_df))
  dummy_0_df <- dummy_0_df |>
    mutate(across(everything(), .fns = identity)) |>
    dplyr::select(all_of(all_cols)) |>
    mutate(across(everything(), ~replace_na(.x, 0)))
  dummy_1_df <- dummy_1_df |>
    mutate(across(everything(), .fns = identity)) |>
    dplyr::select(all_of(all_cols)) |>
    mutate(across(everything(), ~replace_na(.x, 0)))

  X0_cat_encoded[[col]] <- dummy_0_df
  X1_cat_encoded[[col]] <- dummy_1_df
}

# library(proxy)
num_dist <- proxy::dist(x = X0_num, y = X1_num, method = "Euclidean")
num_dist <- as.matrix(num_dist)

cat_dists <- list()
for (col in cat_cols) {
  mat_0 <- as.matrix(X0_cat_encoded[[col]])
  mat_1 <- as.matrix(X1_cat_encoded[[col]])
  dist_mat <- proxy::dist(x = mat_0, y = mat_1, method = "Euclidean")
  cat_dists[[col]] <- as.matrix(dist_mat)
}

combined_cost <- num_dist
for (i in seq_along(cat_dists)) {
  combined_cost <- combined_cost + cat_dists[[i]]
}

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

transport_plan <- matrix(0, nrow = n0, ncol = n1)
for(i in seq_len(nrow(transport_res))) {
  transport_plan[transport_res$from[i], transport_res$to[i]] <- transport_res$mass[i]
}

num_transported <- n0 * (transport_plan %*% as.matrix(X1_num))

cat_transported <- list()
for (col in cat_cols) {
  cat_probs <- transport_plan %*% as.matrix(X1_cat_encoded[[col]])
  cat_encoded_columns <- colnames(X1_cat_encoded[[col]])
  # For each obs., we take the index with the maximum value (approx. proba)
  max_indices <- apply(cat_probs, 1, which.max)
  prefix_pattern <- paste0("^", col, "\\.")
  cat_transported[[col]] <- sapply(
    max_indices, 
    function(x) sub(prefix_pattern, "", cat_encoded_columns[x])
  )
}

tb_ot_transported <- as_tibble(num_transported)
for (col in cat_cols) {
  tb_ot_transported[[col]] <- cat_transported[[col]]
}

save(tb_ot_transported, file = "../output/ot-synthetic.rda")

# Load tb_ot_transported
load("../output/ot-synthetic.rda")
tb_ot_transported <- tb_ot_transported |>
  mutate(X3 = as.factor(X3))
tb_ot_transported <- as.list(tb_ot_transported)

# Compute transport plan
sinkhorn_transport_res <- T4transport::sinkhornD(
  combined_cost, wx = w0, wy = w1, lambda = 0.1
)

sinkhorn_transport_plan <- sinkhorn_transport_res$plan

num_sinkhorn_transported <- n0 * (sinkhorn_transport_plan %*% as.matrix(X1_num))

cat_sinkhorn_transported <- list()
for (col in cat_cols) {
  cat_probs <- sinkhorn_transport_plan %*% as.matrix(X1_cat_encoded[[col]])
  cat_encoded_columns <- colnames(X1_cat_encoded[[col]])
  # For each obs., we take the index with the maximum value (approx. proba)
  max_indices <- apply(cat_probs, 1, which.max)
  prefix_pattern <- paste0("^", col, "\\.")
  cat_sinkhorn_transported[[col]] <- sapply(
    max_indices, 
    function(x) sub(prefix_pattern, "", cat_encoded_columns[x])
  )
}

tb_sinkhorn_transported <- as_tibble(num_sinkhorn_transported)
for (col in cat_cols) {
  tb_sinkhorn_transported[[col]] <- cat_sinkhorn_transported[[col]]
}

save(tb_sinkhorn_transported, file = "../output/sinkhorn-synthetic.rda")

# Load tb_sinkhorn_transported
load("../output/sinkhorn-synthetic.rda")
tb_sinkhorn_transported <- tb_sinkhorn_transported |>
  mutate(X3 = as.factor(X3))
tb_sinkhorn_transported <- as.list(tb_sinkhorn_transported)

library(pbapply)
library(parallel)
ncl <- detectCores()-1
(cl <- makeCluster(ncl))

clusterEvalQ(cl, {
  library(transportsimplex)
  source("../scripts/functions.R")
}) |>
  invisible()

sequential_transport <- seq_trans(
  data = tb_estim, 
  adj = adj, 
  s = A_name, 
  S_0 = 0, # source: untreated
  y = Y_name, 
  num_neighbors = 50, 
  num_neighbors_q = NULL,
  silent = FALSE,
  method = "shortsimplex",
  cl = cl
)

save(sequential_transport, file = "../output/seq-t-synthetic.rda")

stopCluster(cl)

load("../output/seq-t-synthetic.rda")

# library(mediation) # we do not load it
# otherwise it masks a lot of useful functions

# We encode the categorical variable as for optimal transport
tb_med <- tb_estim |> 
  mutate(
    X3 = case_when(
      X3 == "A" ~ 0,
      X3 == "B" ~ 1,
      X3 == "C" ~ 2
    )
  )

med_mod_X1 <- mediation::multimed(
  outcome = "Y", 
  med.main = "X1", 
  med.alt = c("X2", "X3"), 
  treat = "A", 
  data = tb_med
)
# Indirect effect for X1: A -> X1 -> Y
delta_0_med_X1 <- mean((med_mod_X1$d0.lb + med_mod_X1$d0.ub) / 2)
# Direct + Other indirect effects: A -> Y, A -> X2 -> Y, A -> X3 -> Y, 
# A -> X1 -> X2 -> Y, A -> X2 -> X3 -> Y, A -> X1 -> X3 -> Y, A -> X1 -> X2 -> X3 -> Y
zeta_1_med_X1 <- mean((med_mod_X1$z1.lb + med_mod_X1$z1.ub) / 2)
# Total effect
tot_effect_med_X1 <- delta_0_med_X1 + zeta_1_med_X1

med_mod_X2 <- mediation::multimed(
  outcome = "Y", 
  med.main = "X2", 
  med.alt = c("X1", "X3"), 
  treat = "A", 
  data = tb_med
)
# Indirect effect for X2: A -> X2 -> Y, A -> X1 -> X2 -> Y
delta_0_med_X2 <- mean((med_mod_X2$d0.lb + med_mod_X2$d0.ub) / 2)
# Direct + Other indirect effects: A -> Y, A -> X1 -> Y, A -> X3 -> Y, 
# A -> X1 -> X3 -> Y, A -> X2 -> X3 -> Y, A -> X1 -> X2 -> X3
zeta_1_med_X2 <- mean((med_mod_X2$z1.lb + med_mod_X2$z1.ub) / 2)
# Total effect
tot_effect_med_X2 <- delta_0_med_X2 + zeta_1_med_X2

med_mod_X3 <- mediation::multimed(
  outcome = "Y", 
  med.main = "X3", 
  med.alt = c("X1", "X2"), 
  treat = "A", 
  data = tb_med
)
# Indirect effect for ccd: A -> X3 -> Y, A -> X1 -> X3 -> Y, A -> X2 -> X3 -> Y,
# A -> X1 -> X2 -> X3 -> Y
delta_0_med_X3 <- mean((med_mod_X3$d0.lb + med_mod_X3$d0.ub) / 2)
# Direct + Other indirect effects: A -> Y, A -> X1 -> Y, A -> X2 -> Y, 
# A -> X1 -> X2 -> Y
zeta_1_med_X3 <- mean((med_mod_X3$z1.lb + med_mod_X3$z1.ub) / 2)
# Total effect
tot_effect_med_X3 <- delta_0_med_X3 + zeta_1_med_X3

# Total effect
tot_effect_med <- tot_effect_med_X1
# Indirect effects
delta_0_med <- delta_0_med_X1 + delta_0_med_X2 + delta_0_med_X3
# Direct effect 
zeta_1_med <- tot_effect_med - delta_0_med
cbind(delta_0 = delta_0_med, zeta_1 = zeta_1_med, tot_effect = tot_effect_med)

library(randomForest)

causal_effects_ot <- causal_effects_cf(
  data_untreated = tb_untreated, 
  data_treated = tb_treated,
  data_cf_untreated = as_tibble(tb_ot_transported), 
  Y_name = Y_name, 
  A_name = A_name, 
  A_untreated = A_untreated # 0
)

cbind(
  delta_0 = causal_effects_ot$delta_0,
  zeta_1 = causal_effects_ot$zeta_1, 
  tot_effect = causal_effects_ot$tot_effect
)

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

library(tikzDevice)
source("../scripts/utils.R")

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
      str_detect(var_name, "^delta_0") ~ "$\\delta_i(0)$",
      str_detect(var_name, "^zeta_1") ~ "$\\zeta_i(1)$",
      str_detect(var_name, "^tot_effect") ~ "$\\tau_i(1)$",
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
  # "OT" = "#CC79A7",
  "OT-M" = "#009E73",
  "skh" = "darkgray",
  "seq_1" = "#0072B2", 
  "seq_2" = "#D55E00"
)

export_tikz <- FALSE


file_name <- "xp-simulated-indiv-effects"
width_tikz <- 2.3
height_tikz <- 1.7
if (export_tikz == TRUE)
  tikz(paste0("figs/", file_name, ".tex"), width = width_tikz, height = height_tikz)

layout(
  matrix(1:9, byrow = TRUE, ncol = 3),
  widths = c(1, rep(.9, 2)), heights = c(1, rep(.72, 2))
)

for (i in 1:3) {
  x <- case_when(
    i == 1 ~ causal_effects_ot,
    i == 2 ~ causal_effects_sink_ot,
    i == 3 ~ causal_effects_st
  )
  method <- case_when(
    i == 1 ~ "OT-M",
    i == 2 ~ "SKH",
    i == 3 ~ "ST"
  )
  
  for (var_name in c("delta_0_i", "zeta_1_i", "tot_effect")) {
    mar_bottom <- ifelse(i == 3, 2.1, .6)
    mar_left <- ifelse(var_name == "delta_0_i", 3.1, 2.1)
    mar_top <- ifelse(i == 1, 2.1, .1)
    mar_right <- .4
    printed_method <- ifelse(var_name == "delta_0_i", method, "")
    
    par(mar = c(mar_bottom, mar_left, mar_top, mar_right))
    x_lim_list <- list(
      "delta_0_i" = c(-.6, .6),
      "zeta_1_i" = c(-.6, .6),
      "tot_effect" = c(-1, 1)
    )
    
    plot_hist_effects(
      x = x, var_name = var_name, tikz = export_tikz, 
      fill = colour_methods[i],
      printed_method = printed_method, 
      x_lim = x_lim_list[[var_name]],
      print_main = i == 1,
      print_x_axis = i == 3
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

seeds <- 1:200
res_simul_opt_trans <- vector(mode = "list", length = length(seeds))
res_simul_sink_trans <- vector(mode = "list", length = length(seeds))
res_simul_seq_trans <- vector(mode = "list", length = length(seeds))
res_simul_effects <- vector(mode = "list", length = length(seeds))

# # This chunk is not evaluated.
# # The results from previously run simulations are loaded after this chunk.
# # Each of the 200 replications takes about 9 seconds per run on a MB Pro with
# # a Apple M2 Pro chip and 32GB RAM.
# library(pbapply)
# library(parallel)
# ncl <- detectCores()-1
# (cl <- makeCluster(ncl))
# 
# clusterEvalQ(cl, {
#   library(transportsimplex)
#   source("../scripts/functions.R")
# }) |>
#   invisible()
# 
# for (i in 1:length(seeds)) {
#   cat(paste0("Simulation ", i, "/", length(seeds), "\n"))
#   seed <- seeds[i]
#   tb_all <- gen_data(seed)
#   tb <- tb_all |> select(A, X1, X2, X3, Y) |>
#     mutate(across(where(is.character), ~as.factor(.x)))
# 
#   A_name <- "A"
#   A_untreated <- 0
#   Y_name <- "Y"
# 
#   # Optimal Transport
#   tb_ot_transported <- optimal_transport_cf(
#     tb,
#     Y_name,
#     A_name,
#     A_untreated
#   )
#   tb_ot_transported <- tb_ot_transported |>
#     mutate(X3 = as.factor(X3))
#   tb_ot_transported <- as.list(tb_ot_transported)
#   res_simul_opt_trans[[i]] <- tb_ot_transported
# 
#   # Penalized (0.1) Optimal Transport with Sinkhorn
#   tb_sinkhorn_transported <- optimal_transport_cf(
#     tb,
#     Y_name,
#     A_name,
#     A_untreated,
#     pen = 0.1
#   )
#   tb_sinkhorn_transported <- tb_sinkhorn_transported |>
#     mutate(X3 = as.factor(X3))
#   tb_sinkhorn_transported <- as.list(tb_sinkhorn_transported)
#   res_simul_sink_trans[[i]] <- tb_sinkhorn_transported
# 
#   # Sequential Transport
#   sequential_transport <- seq_trans(
#     data = tb,
#     adj = adj,
#     s = A_name,
#     S_0 = 0, # source: untreated
#     y = Y_name,
#     num_neighbors = 50,
#     num_neighbors_q = NULL,
#     silent = FALSE,
#     cl = cl
#   )
#   res_simul_seq_trans[[i]] <- sequential_transport$transported
# 
#   tb_untreated <- tb |> filter(!!sym(A_name) == !!A_untreated)
#   tb_treated <- tb |> filter(!!sym(A_name) != !!A_untreated)
# 
#   n_untreated <- nrow(tb_untreated)
#   n_treated <- nrow(tb_treated)
# 
#   ## Measuring Causal Effect----
# 
#   # With Causal Mediation Analysis
#   tb_med <- tb |>
#     mutate(
#       X3 = case_when(
#         X3 == "A" ~ 0,
#         X3 == "B" ~ 1,
#         X3 == "C" ~ 2
#       )
#     )
# 
#   med_mod_X1 <- mediation::multimed(
#     outcome = "Y",
#     med.main = "X1",
#     med.alt = c("X2", "X3"),
#     treat = "A",
#     data = tb_med
#   )
#   delta_0_med_X1 <- mean((med_mod_X1$d0.lb + med_mod_X1$d0.ub) / 2)
#   zeta_1_med_X1 <- mean((med_mod_X1$z1.lb + med_mod_X1$z1.ub) / 2)
#   tot_effect_med_X1 <- delta_0_med_X1 + zeta_1_med_X1
# 
#   med_mod_X2 <- mediation::multimed(
#     outcome = "Y",
#     med.main = "X2",
#     med.alt = c("X1", "X3"),
#     treat = "A",
#     data = tb_med
#   )
#   delta_0_med_X2 <- mean((med_mod_X2$d0.lb + med_mod_X2$d0.ub) / 2)
#   zeta_1_med_X2 <- mean((med_mod_X2$z1.lb + med_mod_X2$z1.ub) / 2)
#   tot_effect_med_X2 <- delta_0_med_X2 + zeta_1_med_X2
# 
#   med_mod_X3 <- mediation::multimed(
#     outcome = "Y",
#     med.main = "X3",
#     med.alt = c("X1", "X2"),
#     treat = "A",
#     data = tb_med
#   )
#   delta_0_med_X3 <- mean((med_mod_X3$d0.lb + med_mod_X3$d0.ub) / 2)
#   zeta_1_med_X3 <- mean((med_mod_X3$z1.lb + med_mod_X3$z1.ub) / 2)
#   tot_effect_med_X3 <- delta_0_med_X3 + zeta_1_med_X3
#   # Summary of the causal effects
#   tot_effect_med <- tot_effect_med_X1
#   delta_0_med <- delta_0_med_X1 + delta_0_med_X2 + delta_0_med_X3
#   zeta_1_med <- tot_effect_med - delta_0_med
# 
#   # With Counterfactual Values
#   ## Optimal Transport
#   causal_effects_ot <- causal_effects_cf(
#     data_untreated = tb_untreated,
#     data_treated = tb_treated,
#     data_cf_untreated = as_tibble(tb_ot_transported),
#     Y_name = Y_name,
#     A_name = A_name,
#     A_untreated = A_untreated
#   )
#   ## Sinkhorn Optimal Transport
#   causal_effects_sink_ot <- causal_effects_cf(
#     data_untreated = tb_untreated,
#     data_treated = tb_treated,
#     data_cf_untreated = as_tibble(tb_sinkhorn_transported),
#     Y_name = Y_name,
#     A_name = A_name,
#     A_untreated = A_untreated
#   )
#   ## Sequential Transport
#   causal_effects_st <- causal_effects_cf(
#     data_untreated = tb_untreated,
#     data_treated = tb_treated,
#     data_cf_untreated = as_tibble(sequential_transport$transported),
#     Y_name = Y_name,
#     A_name = A_name,
#     A_untreated = A_untreated
#   )
# 
#   res_simul_effects[[i]] <- tibble(
#     seed = seed,
#     delta_0_med = delta_0_med,
#     zeta_1_med = zeta_1_med,
#     tot_effect_med = tot_effect_med,
#     delta_0_ot = causal_effects_ot$delta_0,
#     zeta_1_ot = causal_effects_ot$zeta_1,
#     tot_effect_ot = causal_effects_ot$tot_effect,
#     delta_0_sink_ot = causal_effects_sink_ot$delta_0,
#     zeta_1_sink_ot = causal_effects_sink_ot$zeta_1,
#     tot_effect_sink_ot = causal_effects_sink_ot$tot_effect,
#     delta_0_st = causal_effects_st$delta_0,
#     zeta_1_st = causal_effects_st$zeta_1,
#     tot_effect_st = causal_effects_st$tot_effect
#   )
# }
# 
# save(
#   res_simul_opt_trans,
#   res_simul_sink_trans,
#   res_simul_seq_trans,
#   res_simul_effects,
#   file = "../output/res_simul.rda"
# )
# 
# stopCluster(cl)

load("../output/res_simul.rda")

causal_effects <- list_rbind(res_simul_effects)

library(stringr)
library(tikzDevice)
export_pdf <- FALSE

labels_strip <- c(
  "$\\bar{\\delta}(0)$",
  "$\\bar{\\zeta}(1)$",
  "$\\bar{\\tau}$"
)
if (export_pdf == FALSE) {
  labels_strip = latex2exp::TeX(labels_strip)
}


data_plot <- causal_effects |> 
  pivot_longer(
    cols = -seed, 
    names_to = "name", values_to = "tau"
  ) |> 
  mutate(
    type = case_when(
      str_detect(name, "^delta") ~ "delta",
      str_detect(name, "^zeta") ~ "zeta",
      str_detect(name, "^tot_effect") ~ "tot_effect",
      TRUE ~ NA_character_
    ),
    type = factor(
      type, 
      levels = c("delta", "zeta", "tot_effect"),
      labels = labels_strip
    ),
    Method = case_when(
      str_detect(name, "_med$") ~ "causal_med",
      str_detect(name, "_st$") ~ "seq_ot",
      str_detect(name, "_sink_ot$") ~ "sink_ot",
      str_detect(name, "_ot$") ~ "ot",
      TRUE ~ NA_character_
    ),
    Method = factor(
      Method,
      levels = c(
        "seq_ot", "sink_ot", "ot", "causal_med"
      ),
      labels = c("ST", "SKH", "OT-M", "CM")
    )
  )
p <- ggplot(
  data = data_plot,
  mapping = aes(x = tau, y = Method)
) +
  geom_violin(
    mapping = aes(fill = Method),
    quantiles =  c(.25, .5, .75),
    quantile.linetype = "solid"
    ) +
  labs(x = NULL, y = NULL)

if (export_pdf == TRUE) {
  p <- p + 
    facet_wrap(
      ~ type, scales = "free_x"
    )
} else {
  p <- p +
    facet_wrap(
      ~ type, scales = "free_x",
      labeller = as_labeller(latex2exp::TeX, default = label_parsed)
    )
}

p <- p +
  scale_fill_manual(
    NULL, 
    values = c(
      "CM" = "#56B4E9",
      "OT-M" = colour_methods[["OT-M"]], 
      "SKH" = colour_methods[["skh"]], 
      "ST" =  colour_methods[["seq_1"]]
    ),
    guide = "none"
  ) +
  theme_paper()

p

if (export_pdf == TRUE) {
  ggplot2_to_pdf(
    plot = p + theme(panel.spacing = unit(0.4, "lines")) +
      scale_x_continuous(
        labels = function(x) paste0("$", x, "$"),
        breaks = scales::pretty_breaks(n = 3)
      ),
    filename = "xp-simulated-violin-mc", path = "figs/", 
    width = 1.27*3.3, height = 1.27*1.3,
    crop = TRUE
  )
  
  system(paste0("pdfcrop figs/xp-simulated-violin-mc.pdf figs/xp-simulated-violin-mc.pdf"))
}
