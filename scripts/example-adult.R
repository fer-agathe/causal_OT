library(tidyverse)
library(seqtransfairness)
library(transportsimplex)
library(cluster)

#' @importFrom cluster daisy
#' @importFrom Hmisc wtd.quantile
#' @importFrom nnet multinom
#' @importFrom purrr map_chr

source("functions.R")
vars <- c(
  "sex", "age", "native_country", "marital_status", "education_num",
  "workclass", "hours_per_week", "occupation", "income"
)
s <- "sex"
y <- "income"

# reading in the UCI Adult data from {fairadapt}
adult <- readRDS(
  system.file("extdata", "uci_adult.rds", package = "fairadapt")
) |>
  as_tibble() |> 
  select(!!vars) |> 
  mutate(
   occupation = fct_lump_prop(occupation, prop = 0.02)
  )

# Test
adult <- adult |> 
  filter(
    !(sex == "Male" & occupation == "Sales")
  ) |> 
  filter(
    !(sex == "Female" & occupation == "Craft-repair")
  )

# Adjacency matrix
adj_mat <- c(
  0, 0, 0, 1, 1, 1, 1, 1, 1, # sex
  0, 0, 0, 1, 1, 1, 1, 1, 1, # age
  0, 0, 0, 1, 1, 1, 1, 1, 1, # native_country
  0, 0, 0, 0, 1, 1, 1, 1, 1, # marital_status
  0, 0, 0, 0, 0, 1, 1, 1, 1, # education_num
  0, 0, 0, 0, 0, 0, 0, 0, 1, # workclass
  0, 0, 0, 0, 0, 0, 0, 0, 1, # hours_per_week
  0, 0, 0, 0, 0, 0, 0, 0, 1, # occupation
  0, 0, 0, 0, 0, 0, 0, 0, 0  # income
) |> matrix(
  nrow = length(vars), ncol = length(vars),
  dimnames = list(vars, vars), byrow = TRUE
)

causal_graph <- fairadapt::graphModel(adj_mat)
plot(causal_graph)

seed <- 2025
# 
# data = adult
# adj = adj_mat
# s = "sex"
# S_0 = "Female"
# y = "income"
# num_neighbors = 50
# num_neighbors_q = NULL
# silent = FALSE

sequential_transport <- seq_trans(
  data = adult, 
  adj = adj_mat, 
  s = "sex", 
  S_0 = "Female", 
  y = "income", 
  num_neighbors = 50, 
  num_neighbors_q = NULL,
  silent = FALSE
)

save(sequential_transport, file = "sequential_transport_adult.rda")


