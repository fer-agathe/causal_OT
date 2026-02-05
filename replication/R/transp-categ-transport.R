# col_categ <- c("#ffdd55","#944edf","#3fb3b2")
col_group <- c("#00A08A","#F2AD00", "#1b95e0")
col_categ <- c("#56B4E9", "#D55E00", "#CC79A7")
colA <- col_categ[1] ; colB <- col_categ[2] ; colC <- col_categ[3]
colGpe1 <- col_group[2]
colGpe0 <- col_group[1]
font_size <- 20
font_family <- "CMU Serif"

path <- "./figs/"

theme_ggtern_paper <- function(...) {
  font_family <- "CMU Serif"
  font_size <- 10
  theme(
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text.x = element_text(colour = "black"),
    strip.text = ggtext::element_markdown(),
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    axis.title = element_text(size = rel(1)),
    tern.axis.arrow.show = TRUE,
    tern.axis.arrow.sep = .13,
    tern.axis.vshift = .05,
    panel.border = element_rect(colour = NA)
  )
}

theme_ggtern_minimal <- function(...) {
  font_family <- "CMU Serif"
  font_size <- 10
  theme(
    strip.background = element_rect(colour = "black", fill = NA),
    # strip.text.x = element_text(colour = "black"),
    # strip.text = ggtext::element_markdown(),
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    axis.title = element_text(size = rel(1)),
    tern.axis.arrow.show = FALSE,
    tern.axis.arrow.sep = .13,
    tern.axis.vshift = .05,
    panel.border = element_rect(colour = NA),
    tern.axis.text.T = element_blank(),
    tern.axis.text.L = element_blank(),
    tern.axis.text.R = element_blank(),
    tern.axis.vshift = 0.2
  )
}


#' Theme for ggplot2
#'
#' @param ... Arguments passed to the theme function.
#' @export
#' @importFrom ggplot2 element_rect element_text element_blank element_line unit
#'   rel theme
#'
theme_paper <- function (...) {
  theme(
    text = element_text(family = font_family),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    # panel.border = element_rect(fill = NA, colour = "black", linewidth = 1),
    panel.border = element_blank(),
    # axis.line = element_line(color = "black"),
    axis.line <- element_blank(),
    axis.text = element_text(color = "black"),
    legend.text = element_text(size = rel(1)),
    legend.title = element_text(size = rel(1)),
    legend.background = element_rect(fill = "transparent", color = NULL),
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    # legend.box = "vertical",
    legend.key = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_line(colour = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0, size = rel(1), face = "bold"),
    plot.title.position = "plot",
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text = element_text(size = rel(1))
  )
}

library(dplyr)
library(tidyr)
library(stringr)
library(ggtern)
library(compositions)

set.seed(1234)
n <- 100
n0 <- 100
n1 <- 110
p0 <- c(0.1, 0.5, 0.4)
p1 <- c(0.5, 0.3, 0.2)
# Sample category
x0 <- sample(c(rep("A", p0[1] * n0), rep("B", p0[2] * n0), rep("C", p0[3] * n0)), replace = FALSE)
x1 <- sample(c(rep("A", p1[1] * n1), rep("B", p1[2] * n1), rep("C", p1[3] * n1)), replace = FALSE)
cat_levels <- c("A", "B", "C")

source("../scripts/utils.R")

p_barplot_source <- ggplot(
  data = tibble(x0 = x0) |> count(x0) |> 
    arrange(desc(x0)) |> 
    mutate(
      prop = n / sum(n),
      lab_y = cumsum(prop) - prop/2
    )
) +
  geom_bar(stat = "identity", mapping = aes(x = factor(1), y = prop, fill = x0)) +
  geom_text(mapping = aes(x = factor(1), y = lab_y, label = x0), family = font_family) +
  scale_fill_manual(values = col_categ, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_blank(), 
    panel.grid.major.x = element_blank()
  )
p_barplot_source

library(MCMCpack)
set.seed(12345)
alpha_A <- c(9, 3, 2)
Z_A <- as.data.frame(rdirichlet(n0 + n1, alpha_A))
alpha_B <- c(3, 11, 4)
Z_B <- as.data.frame(rdirichlet(n0 + n1, alpha_B))
alpha_C <- c(2, 3, 9)
Z_C <- as.data.frame(rdirichlet(n0 + n1, alpha_C))
# For each observation from group 0 and matched obs from group 1, we have
# drawn a category (A, B, or C).
# We add drawn propensities, depending on the category
Z <- Z_A
category <- c(x0, x1)
Z[category == "B", ] <- Z_B[category == "B", ]
Z[category == "C", ] <- Z_C[category == "C", ]
tb_sample_z <- as_tibble(Z)
names(tb_sample_z) <- c("A", "B", "C")
tb_sample_z$group <- factor(c(rep(0, n0), rep(1, n1)), levels = c(0, 1))

tb_sample_z_1 <- tb_sample_z

p_compositional <- ggtern(
  data = tb_sample_z_1 |> filter(group == 0),
  mapping = aes(x = A, y = B, z = C)) +
  geom_point(alpha = .8, size = .5, mapping = aes(color = group)) +
  scale_colour_manual(name = "group",values = col_group, guide = "none") +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_minimal()
  
p_compositional

#' Pairwise distance matrix on the simplex
#'
#' @description
#' Computes the pairwise distance matrix of observations in the simplex, using
#' the cost function for optimal transport on the unit simplex as the distance
#' metric.
#'
#' @param X Matrix of observations (one observation per row).
#' @param Y Matrix of observations (one observation per row).
#'
#' @returns A matrix of size n x m, where n is the number of observation in X,
#'  and m is the number of observations in Y, containing the distances between
#'  observations in X and Y.
#' @noRd
compute_pdist_simplex_fast <- function(X, Y) {
  
  p <- ncol(X)
  invX <- 1 / X
  
  # R[j,i] = sum_k Y[j,k] * invX[i,k]
  R <- Y %*% t(invX)
  
  logXmean <- rowMeans(log(X))
  logYmean <- rowMeans(log(Y))
  
  # M[i,j] = log(R[j,i]) - log(p) - logYmean[j] + logXmean[i]
  M_t <- log(R) - log(p) -
    outer(logYmean, rep(1, length(logXmean))) +
    outer(rep(1, length(logYmean)), logXmean)
  
  t(M_t)
}

#' Solving the Optimal Transport Problem
#'
#' @description
#' Finds the optimal transport plan using shortsimplex method.
#'
#' @param dxy Cost matrix of transport distances between points in X and Y.
#' @param wx Weights (marginal distribution) for X.
#' @param wy Weights (marginal distribution) for Y.
#' @param p Order of the Wassterstein distance. (If p=2: squared Euclidean
#'  cost).
#'
#' @importFrom transport transport
#'
#' @noRd
wass_lp_fast <- function(dxy, 
                         wx, 
                         wy, 
                         p = 2) {
  
  stopifnot(all(abs(sum(wx) - 1) < 1e-8), all(abs(sum(wy) - 1) < 1e-8))
  
  m <- length(wx)
  n <- length(wy)
  
  # Convert dxy to a cost matrix (flattened)
  cost <- as.matrix(dxy)^p
  
  # Solve the OT problem (default method = "shortsimplex")
  plan <- transport::transport(wx, wy, costm = cost)
  
  # Convert transport plan (sparse format) to matrix
  gamma <- matrix(0, m, n)
  for (i in seq_len(nrow(plan))) {
    gamma[plan$from[i], plan$to[i]] <- plan$mass[i]
  }
  
  # Compute Wasserstein distance
  value <- sum(gamma * cost)^(1 / p)
  
  list(distance = value, plan = gamma)
}

tb_sample_z_1_0 <- tb_sample_z_1 |> filter(group == 0)
tb_sample_z_1_1 <- tb_sample_z_1 |> filter(group == 1)

prop_0 <- tb_sample_z_1_0 |> 
  dplyr::select(-group) |> 
  as.matrix()

prop_1 <- tb_sample_z_1_1 |> 
  dplyr::select(-group) |> 
  as.matrix()

dist_mat <- compute_pdist_simplex_fast(X = prop_0, Y = prop_1)

par_w0 <- rep(1/nrow(prop_0), nrow(prop_0))
par_w1 <- rep(1/nrow(prop_1), nrow(prop_1))
# Solve the optimal transport problem
ot_plan <- wass_lp_fast(dxy = dist_mat, wx = par_w0, wy = par_w1, p = 2)

ot_plan_tb <- tibble(
  from = 1:n0,
  to = max.col(ot_plan$plan, ties.method = "first")
)

p_matching_1 <- ggtern(
  data = tb_sample_z_1,
  mapping = aes(x = A, y = B, z = C)) +
  geom_point(size = .5, alpha = 0.8, mapping = aes(color = group)) +
  scale_colour_manual(name = "group",values = col_group, guide = "none") +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_minimal()


# Create interpolated values using McCann (1997) displacement
f_line_simplex <- function(x, 
                           y, 
                           lgt = 601) {
  
  zx <- as.numeric(clr(x))[1:2]
  zy <- as.numeric(clr(y))[1:2]
  t <- seq(0, 1, length = lgt)
  
  tx <- cbind(
    (1 - t) * zx[1] + t * zy[1], 
    (1 - t) * zx[2] + t * zy[2]
  )
  tx <- cbind(tx, -(tx[, 1] + tx[, 2]))
  df <- as.data.frame(matrix(as.numeric(clrInv(tx)), lgt, 3))
  names(df) <- c("A","B","C")
  
  df
}


for (i in 1:nrow(ot_plan_tb)) {
  lines_1 <- f_line_simplex(
    x = prop_0[i, 1:3], 
    y = prop_1[ot_plan_tb$to[i], 1:3], 
    lgt = 101
  ) |> 
    as_tibble()
  
  p_matching_1 <- p_matching_1 + 
    geom_line(
      data = lines_1, 
      mapping = aes(x = A, y = B, z = C), 
      color = col_group[1], linewidth = .2,, alpha = .5,
      arrow = arrow(length = unit(0.10, "cm"))
    )
}

p_matching_1

p_ggtern_matched_1 <- ggtern(
  data = tb_sample_z_1 |> filter(group == 1),
  mapping = aes(x = A, y = B, z = C)
) +
  geom_point(alpha = .8, size = .5, colour = col_group[2]) +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_minimal()
  
p_ggtern_matched_1

set.seed(123)
n <- 1000
# Concentration parameter:
alpha <-  c(2, 5, 3)
# Draw n samples
samples <- rdirichlet(n, alpha = alpha)

# Unit vectors of S_3
vertices <- matrix(c(
  1, 0, 0,  # A
  0, 1, 0,  # B
  0, 0, 1   # C
), byrow = TRUE, ncol = 3)

# source weights
mass_source <- rep(1 / n, n)
# target weights
mass_target <- c(3, 2, 1) / 6

# Cost matrix (squared Euclidean distance)
cost_matrix <- as.matrix(dist(rbind(samples, vertices))^2)
cost_matrix <- cost_matrix[1:n, (n + 1):(n + 3)]

# We assign eah observation to one vertex
# by minimizing the global transport cost, while matching marginals

# Solve the optimal transport plan
ot_plan <- transport::transport(
  a = mass_source, b = mass_target, costm = cost_matrix, 
  method = "shortsimplex"
)

# Assign each sample to a category based on OT plan
assignment <- rep(NA, n)
# mass each source sends to each target
mass_matrix <- matrix(0, nrow = n, ncol = 3)

for (j in 1:nrow(ot_plan)) {
  from <- ot_plan$from[j]
  to <- ot_plan$to[j]
  mass <- ot_plan$mass[j]
  mass_matrix[from, to] <- mass_matrix[from, to] + mass
}

# Assign each source point to the target it contributes the most mass to
assignment <- max.col(mass_matrix, ties.method = "first")

colnames(samples) <- c("A", "B", "C")

samples_dirichlet_p <- 
  as_tibble(samples) |> 
  mutate(category = colnames(samples)[assignment])
samples_dirichlet_p

ggtern(
  data = samples_dirichlet_p,
  mapping = aes(x = A, y = B, z = C, color = category)
) +
  geom_point(alpha = .8, size = .5) +
  scale_colour_manual(name = "category", values = col_categ) +
  theme_light(base_size = font_size, base_family = font_family) +
  theme(
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text.x = element_text(colour = "black"),
    strip.text.y = element_text(colour = "black"),
    # strip.text = ggtext::element_markdown(),
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    axis.title = element_text(size = rel(.8)),
    tern.axis.arrow.show = TRUE,
    tern.axis.arrow.sep = .16,
    tern.axis.vshift = .09,
    legend.position = "bottom",
    legend.title = element_text(size = .8 * font_size),
    legend.text = element_text(size = .8 * font_size),
    panel.border = element_rect(colour = NA)
  ) +
  theme_hidetitles() +
  guides(colour = guide_legend(override.aes = list(size = 2)))


generate_simplex_grid <- function(resolution = 100) {
  grid <- expand_grid(
    A = seq(0, 1, length.out = resolution),
    B = seq(0, 1, length.out = resolution)
  ) |> 
    mutate(C = 1 - A - B) |> 
    filter(C >= 0)
  
  as_tibble(grid)
}

get_category_density_2D <- function(samples, 
                                    grid_points,
                                    category_label) {
  
  category_data <- samples |> 
    filter(category == !!category_label) |> 
    dplyr::select("A", "B")
  eval_grid <- grid_points[, c("A", "B")]
  
  if (nrow(category_data) < 10) {
    H <- diag(0.01, 2)
  } else {
    H <- tryCatch(Hpi(category_data), error = function(e) diag(0.01, 2))
  }
  
  kde_result <- ks::kde(x = category_data, H = H, eval.points = eval_grid)
  
  kde_result$estimate
}

grid_points <- generate_simplex_grid(resolution = 100)

dens_A <- get_category_density_2D(
  samples = samples_dirichlet_p, grid_points = grid_points, category_label = "A"
)
dens_B <- get_category_density_2D(
  samples = samples_dirichlet_p, grid_points = grid_points, category_label = "B"
)
dens_C <- get_category_density_2D(
  samples = samples_dirichlet_p, grid_points = grid_points, category_label = "C"
)

min_dens <- pmin(dens_A, dens_B, dens_C)
max_idx <- which.max(min_dens)
intersection_point <- grid_points[max_idx, ]
tb_intersection <- as_tibble(intersection_point)

p <- ggtern(
    data = samples_dirichlet_p,
  mapping = aes(x = A, y = B, z = C)
) +
  geom_point(alpha = .8, size = .5, mapping = aes(color = category)) +
  geom_point(data = tb_intersection) +
  scale_colour_manual(values = col_categ) +
  theme_light(base_size = font_size, base_family = font_family) +
  theme(
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text.x = element_text(colour = "black"),
    strip.text.y = element_text(colour = "black"),
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    axis.title = element_text(size = rel(.8)),
    tern.axis.arrow.show = TRUE,
    tern.axis.arrow.sep = .16,
    tern.axis.vshift = .09,
    legend.position = "bottom",
    legend.title = element_text(size = .8 * font_size),
    legend.text = element_text(size = .8 * font_size),
    panel.border = element_rect(colour = NA)
  ) +
  theme_hidetitles() +
  guides(colour = guide_legend(override.aes = list(size = 2)))

p

#' @param n Number of observations to sample from the Dirichlet Distribution.
#' @param Vector of shape parameters, or matrix of shape parameters 
#'  corresponding to the number of draw. Default to \eqn{(1,1,1)}.
#' @param p Vector of target probabilities. Default to \eqn{(1/3, 1/3, 1/3)}.
#' 
get_data_assignment <- function(n,
                                alpha = c(1, 1, 1),
                                p = c(1, 1, 1) / 3,
                                intersection_point = TRUE) {
  
  # Draw n samples
  samples <- rdirichlet(n, alpha = alpha)
  
  # Unit vectors of S_3
  vertices <- matrix(c(
    1, 0, 0,  # A
    0, 1, 0,  # B
    0, 0, 1   # C
  ), byrow = TRUE, ncol = 3)
  
  # source weights
  mass_source <- rep(1 / n, n)
  # target weights
  mass_target <- p
  
  # Cost matrix (squared Euclidean distance)
  cost_matrix <- as.matrix(dist(rbind(samples, vertices))^2)
  cost_matrix <- cost_matrix[1:n, (n + 1):(n + 3)]
  
  # We assign eah observation to one vertex
  # by minimizing the global transport cost, while matching marginals
  
  # Solve the optimal transport plan
  ot_plan <- transport::transport(
    a = mass_source, b = mass_target, costm = cost_matrix, 
    method = "shortsimplex"
  )
  
  # Assign each sample to a category based on OT plan
  assignment <- rep(NA, n)
  # mass each source sends to each target
  mass_matrix <- matrix(0, nrow = n, ncol = 3)
  
  for (j in 1:nrow(ot_plan)) {
    from <- ot_plan$from[j]
    to <- ot_plan$to[j]
    mass <- ot_plan$mass[j]
    mass_matrix[from, to] <- mass_matrix[from, to] + mass
  }
  
  # Assign each source point to the target it contributes the most mass to
  assignment <- max.col(mass_matrix, ties.method = "first")
  
  colnames(samples) <- c("A", "B", "C")
  samples <- 
    as_tibble(samples) |> 
    mutate(category = colnames(samples)[assignment])
  
  #
  # Intersection point
  #
  if (intersection_point == TRUE) {
    # Create a triangular grid over (A,B,C) constrained to S_3
    grid_points <- generate_simplex_grid(resolution = 100)
    
    # Evaluation of KDE for data from the source distribution
    dens_A <- get_category_density_2D(
      samples = samples, grid_points = grid_points, category_label = "A"
    )
    dens_B <- get_category_density_2D(
      samples = samples, grid_points = grid_points, category_label = "B"
    )
    dens_C <- get_category_density_2D(
      samples = samples, grid_points = grid_points, category_label = "C"
    )
    # Find point where min(densities) is maximal
    min_dens <- pmin(dens_A, dens_B, dens_C)
    max_idx <- which.max(min_dens)
    
    intersection_point <- grid_points[max_idx, ]
    
    tb_intersection <- as_tibble(intersection_point)
  } else {
    tb_intersection <- NULL
  }
  
  list(
    samples = samples,
    tb_intersection = tb_intersection
  )
}

samples_unif_unif <- get_data_assignment(
  n = n, 
  alpha = c(1, 1, 1), 
  p = c(1, 1, 1) / 3
)

samples_unif_p <- get_data_assignment(
  n = n, 
  alpha = c(1, 1, 1), 
  p = c(3, 2, 1) / 6
)

samples_dirichlet_unif <- get_data_assignment(
  n = n, 
  alpha = c(2, 5, 3), 
  p = c(1, 1, 1) / 3
)

samples_dirichlet_p <- get_data_assignment(
  n = n, 
  alpha = c(2, 5, 3), 
  p = c(3, 2, 1) / 6
)

p <- ggtern(
  data = samples_unif_unif$samples |> 
    mutate(distrib = "Uniform", p = "(1,1,1)/3") |> 
    bind_rows(
      samples_unif_p$samples |> 
        mutate(distrib = "Uniform", p = "(3,2,1)/6")
    ) |> 
    bind_rows(
      samples_dirichlet_unif$samples |> 
        mutate(distrib = "2_5_3", p = "(1,1,1)/3")
    ) |> 
    bind_rows(
      samples_dirichlet_p$samples |> 
        mutate(distrib = "2_5_3", p = "(3,2,1)/6")
    ) |>
    mutate(
      distrib = factor(
        distrib,
        levels = c("Uniform", "2_5_3"),
        labels = c(
          "Uniform" = parse(text = latex2exp::TeX("D(1,1,1)")),
          "2_5_3" = parse(text = latex2exp::TeX("$D(2,5,3)$"))
        )
      ),
      p = factor(
        p,
        levels = c("(1,1,1)/3", "(3,2,1)/6"),
        labels = c(
          "(1,1,1)/3" = parse(text = latex2exp::TeX("$p=(1,1,1)/3$")),
          "(3,2,1)/6" = parse(text = latex2exp::TeX("$p=(3,2,1)/6$"))
        )
      )
    ),
  mapping = aes(x = A, y = B, z = C)
) +
  geom_point(alpha = .8, size = .5, mapping = aes(color = category)) +
  # geom_point(
  #   data = samples_unif_unif$tb_intersection |> 
  #     mutate(distrib = "Uniform", p = "(1,1,1)/3") |> 
  #     bind_rows(
  #       samples_unif_p$tb_intersection |> 
  #         mutate(distrib = "Uniform", p = "(3,2,1)/6")
  #     ) |> 
  #     bind_rows(
  #       samples_dirichlet_unif$tb_intersection |> 
  #         mutate(distrib = "2_5_3", p = "(1,1,1)/3")
  #     ) |> 
  #     bind_rows(
  #       samples_dirichlet_p$tb_intersection |> 
  #         mutate(distrib = "2_5_3", p = "(3,2,1)/6")
  #     ) |>
  #     mutate(
  #       distrib = factor(
  #         distrib,
  #         levels = c("Uniform", "2_5_3"),
  #         labels = c(
  #           "Uniform" = parse(text = latex2exp::TeX("D(1,1,1)")),
  #           "2_5_3" = parse(text = latex2exp::TeX("$D(2,5,3)$"))
  #         )
  #       ),
  #       p = factor(
  #         p,
  #         levels = c("(1,1,1)/3", "(3,2,1)/6"),
  #         labels = c(
  #           "(1,1,1)/3" = parse(text = latex2exp::TeX("$p=(1,1,1)/3$")),
  #           "(3,2,1)/6" = parse(text = latex2exp::TeX("$p=(3,2,1)/6$"))
  #         )
  #       )
  #     )
  # ) +
  scale_colour_manual(values = col_categ) +
  facet_grid(p ~ distrib, labeller = label_parsed, switch = "y") +
  theme_light(base_size = font_size, base_family = font_family) +
  theme(
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text.x = element_text(colour = "black"),
    strip.text.y = element_text(colour = "black"),
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    axis.title = element_text(size = rel(.8)),
    tern.axis.arrow.show = TRUE,
    tern.axis.arrow.sep = .16,
    tern.axis.vshift = .09,
    legend.position = "bottom",
    legend.title = element_text(size = .8 * font_size),
    legend.text = element_text(size = .8 * font_size),
    panel.border = element_rect(colour = NA)
  ) +
  theme_hidetitles() +
  guides(colour = guide_legend(override.aes = list(size = 2)))

p

# filename <- "baryc-centr-bal"
# ggsave(
#   p, file = str_c(path, filename, ".pdf"),
#   height = 3.3*1.75, width = 3.25*1.75,
#   family = font_family,
#   device = cairo_pdf
# )
# # Crop PDF
# system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))

#' OT for categorical variable, from source distribution to target 
#' probabilities.
#' 
#' @param probs Propensities from the source distribution (individuals in rows,
#'  classes in columns).
#' @param labels Levels (labels) of the classes.
#' @param p Vector of target probabilities. If omitted, uniform weights are 
#'  used.
#' 
get_assignment <- function(probs,
                           labels,
                           p = NULL) {
  
  n_labels <- ncol(probs)
  n <- nrow(probs)
  if (is.null(p)) p <- rep(1, n_labels) / n_labels # Uniform weights
  
  # Unit vectors
  vertices <- diag(n_labels)
  # colnames(vertices) <- colnames()
  # source weights
  mass_source <- rep(1 / n, n)
  # target weights
  mass_target <- as.numeric(p)
  
  # Cost matrix (squared Euclidean distance)
  cost_matrix <- as.matrix(dist(rbind(probs, vertices))^2)
  cost_matrix <- cost_matrix[1:n, (n + 1):(n + n_labels)]
  
  # Assign each observation to one vertex
  # by minimizing the global transport cost, while matching marginals
  
  # Solve the optimal transport plan
  ot_plan <- transport::transport(
    a = mass_source, b = mass_target, costm = cost_matrix, 
    method = "shortsimplex"
  )
  
  # Assign each sample to a category based on OT plan
  assignment <- rep(NA, n)
  # mass each source sends to each target
  mass_matrix <- matrix(0, nrow = n, ncol = n_labels)
  
  for (j in 1:nrow(ot_plan)) {
    from <- ot_plan$from[j]
    to <- ot_plan$to[j]
    mass <- ot_plan$mass[j]
    mass_matrix[from, to] <- mass_matrix[from, to] + mass
  }
  
  # Assign each source point to the target it contributes the most mass to
  assignments <- max.col(mass_matrix, ties.method = "random")
  #factor(c(1, 2, 4), levels = 1:4, labels = c("A", "B", "C", "D"))
  
  factor(assignments, levels = 1:length(labels), labels = labels)
}

transported <- get_assignment(
  probs = tb_sample_z_1 |> filter(group == 0) |> dplyr::select(-group) |> as.matrix(), 
  labels = c("A", "B", "C"), 
  p = p1
)
head(transported)

ggtern(
  data = tb_sample_z_1 |> filter(group == 0) |> 
    mutate(category = x0) |> 
    bind_rows(
      tb_sample_z_1 |> filter(group == 1) |> 
        mutate(category = x1)
    ) |> 
    mutate(
      group = factor(
        group, 
        levels = c(0, 1), 
        labels = c(
          paste0("<span style='color:", colGpe0,";'>Group 0</span>"), 
          paste0("<span style='color:", colGpe1,";'>Group 1</span>")
        )
      )
    ),
  mapping = aes(x = A, y = B, z = C)
) +
  geom_point(alpha = .8, size = .5, mapping = aes(color = category)) +
  facet_wrap(~ group) +
  scale_colour_manual(values = col_categ, guide = "none") +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_minimal() +
  theme(
    strip.text.x = element_text(colour = "black"),
    strip.text = ggtext::element_markdown()
  )


p_ggtern_assignment_1 <- ggtern(
  data = tb_sample_z_1 |> filter(group == 0) |> 
    mutate(category = transported),
  mapping = aes(x = A, y = B, z = C)
) +
  geom_point(alpha = .8, size = .5, mapping = aes(color = category)) +
  scale_colour_manual(values = col_categ, guide = "none") +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_minimal()
p_ggtern_assignment_1

p_barplot_target_1 <- ggplot(
  data = tibble(x_t = transported) |> count(x_t) |> 
    arrange(desc(x_t)) |> 
    mutate(
      prop = n / sum(n),
      lab_y = cumsum(prop) - prop/2
    )
) +
  geom_bar(stat = "identity", mapping = aes(x = factor(1), y = prop, fill = x_t)) +
  geom_text(mapping = aes(x = factor(1), y = lab_y, label = x_t), family = font_family) +
  scale_fill_manual(values = col_categ, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(
    axis.ticks.x = element_blank(), 
    axis.text.x = element_blank(), 
    panel.grid.major.x = element_blank()
  )
p_barplot_target_1

# export_graph <- FALSE
# 
# if (export_graph) {
# 
#   library(tikzDevice)
#   path <- "figs/"
#   filename <- "transp-categ-barplot-source"
#   ggplot2_to_pdf(
#     plot = p_barplot_source +
#       scale_y_continuous(labels = function(x) paste0("$", x, "$")),
#     filename = filename, path = path,
#     width = .8, height = 1.6,
#     crop = TRUE
#   )
#   system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
# 
#   filename <- "transp-categ-barplot-target"
#   ggplot2_to_pdf(
#     plot = p_barplot_target_1 +
#       scale_y_continuous(labels = function(x) paste0("$", x, "$")),
#     filename = filename, path = path,
#     width = .8, height = 1.6,
#     crop = TRUE
#   )
#   system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
# 
#   filename <- "transp-categ-comp"
#   ggplot2_to_pdf(
#     plot = p_compositional + theme(tern.axis.title.show = FALSE),
#     filename = filename, path = "figs/",
#     width = 3, height = 1.6,
#     crop = TRUE
#   )
#   system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
# 
#   filename <- "transp-categ-matching"
#   ggplot2_to_pdf(
#     plot = p_matching_1 + theme(tern.axis.title.show = FALSE),
#     filename = filename, path = path,
#     width = 3, height = 1.6,
#     crop = TRUE
#   )
#   system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
# 
#   filename <- "transp-categ-matched"
#   ggplot2_to_pdf(
#     plot = p_ggtern_matched_1 + theme(tern.axis.title.show = FALSE),
#     filename = filename, path = path,
#     width = 3, height = 1.6,
#     crop = TRUE
#   )
#   system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
# 
#   filename <- "transp-categ-assignment"
#   ggplot2_to_pdf(
#     plot = p_ggtern_assignment_1 + theme(tern.axis.title.show = FALSE),
#     filename = filename, path = path,
#     width = 3, height = 1.6,
#     crop = TRUE
#   )
#   system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
# }
