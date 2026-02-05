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

# Also required:
# install.packages(mlr3fairness)



library(extrafont, quietly = TRUE)
font_family <- "CMU Serif"


path <- "./figs/"
if (!dir.exists(path)) dir.create(path)

source("../scripts/functions.R")

data(compas, package = "mlr3fairness")

tb <- compas |> 
  as_tibble() |> 
  filter(race %in% c("Caucasian", "African-American")) |>
  select(
    race, # sensitive
    age, 
    priors_count, # The prior criminal records of defendants. 
    c_charge_degree, # F: Felony M: Misdemeanor
    is_recid # outcome
  ) |> 
  mutate(
   race = ifelse(race == "African-American", 0, 1), # African-American as "untreated"
   is_recid = ifelse(is_recid == 0, 0, 1)
  )

dim(tb)
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

seed <- 1234
set.seed(seed)

A_name <- "race" # treatment name
Y_name <- "is_recid" # outcome name
A_untreated <- 0
A <- tb[[A_name]]
ind_untreated <- which(A == A_untreated)
tb_untreated <- tb[ind_untreated, ]
tb_treated <- tb[-ind_untreated, ]

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
  cat_transported[[col]] <- sapply(max_indices, function(x) sub(prefix_pattern, "", cat_encoded_columns[x]))
}

tb_ot_transported <- as_tibble(num_transported)
for (col in cat_cols) {
  tb_ot_transported[[col]] <- cat_transported[[col]]
}

save(tb_ot_transported, file = "../output/ot-compas.rda")

# Load tb_ot_transported
load("../output/ot-compas.rda")
tb_ot_transported <- tb_ot_transported |>
  mutate(c_charge_degree = as.factor(c_charge_degree))
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
  cat_sinkhorn_transported[[col]] <- sapply(max_indices, function(x) sub(prefix_pattern, "", cat_encoded_columns[x]))
}

tb_sinkhorn_transported <- as_tibble(num_sinkhorn_transported)
for (col in cat_cols) {
  tb_sinkhorn_transported[[col]] <- cat_sinkhorn_transported[[col]]
}

save(tb_sinkhorn_transported, file = "../output/sinkhorn-compas.rda")

# Load tb_sinkhorn_transported
load("../output/sinkhorn-compas.rda")
tb_sinkhorn_transported <- tb_sinkhorn_transported |>
  mutate(c_charge_degree = as.factor(c_charge_degree))
tb_sinkhorn_transported <- as.list(tb_sinkhorn_transported)

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

save(sequential_transport, file = "../output/seq-t-compas.rda")

stopCluster(cl)
# Load sequential_transport
load("../output/seq-t-compas.rda")

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

# Total effect
tot_effect_med <- tot_effect_med_age
# Indirect effects
delta_0_med <- delta_0_med_age + delta_0_med_pc + delta_0_med_ccd
# Direct effect 
zeta_1_med <- tot_effect_med - delta_0_med
cbind(delta_0 = delta_0_med, zeta_1 = zeta_1_med, tot_effect = tot_effect_med)

library(randomForest)

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


file_name <- "compas-dist-indiv-effects"
width_tikz <- 2.25
height_tikz <- 1.4
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

tibble::tribble(
  ~method, ~delta_0, ~zeta_1, ~tot_effect,
  "CM", delta_0_med, zeta_1_med, tot_effect_med,
  "OT", causal_effects_ot$delta_0, causal_effects_ot$zeta_1, causal_effects_ot$tot_effect,
  "SKH", causal_effects_sink_ot$delta_0, causal_effects_sink_ot$zeta_1, causal_effects_sink_ot$tot_effect,
  "ST", causal_effects_st$delta_0, causal_effects_st$zeta_1, causal_effects_st$tot_effect
) |> 
  knitr::kable(digits = 2)

# Load sequential_transport
load("../output/seq-t-compas.rda")

causal_effects_st$delta_0

# Observations
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

modified_causal_effects_cf(data_untreated, data_cf_untreated, mu_untreated_model)$delta_0

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
