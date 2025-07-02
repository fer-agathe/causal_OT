# Categorical Representation of Compositional
# Features, or Barycentric Centroid of Balance

library(tidyverse)
library(transport)
library(MCMCpack)
library(ggtern)

set.seed(123)
n <- 1000

library(extrafont)
loadfonts(device = "pdf")
colours <- c("A" = "#FFE788", "B" = "#9F78C4", "C" = "#5FC7C6")
font_size <- 20
font_family <- "CMU Serif"

path <- "./figs/"
if (!dir.exists(path)) dir.create(path)

#' @param n Number of observations to sample from the Dirichlet Distribution.
#' @param Vector of shape parameters, or matrix of shape parameters 
#'  corresponding to the number of draw. Default to \eqn{(1,1,1)}.
#' @param p Vector of target probabilities. Default to \eqn{(1/3, 1/3, 1/3)}.
#' 
get_data_assignment <- function(n,
                                alpha = c(1, 1, 1),
                                p = c(1, 1, 1) / 3) {
  
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
  
  as_tibble(samples) |> 
    mutate(category = colnames(samples)[assignment])
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



p <- 
  ggtern(
  data = samples_unif_unif |> mutate(distrib = "Uniform", p = "(1,1,1)/3") |> 
    bind_rows(
      samples_unif_p |> 
        mutate(distrib = "Uniform", p = "(3,2,1)/6")
    ) |> 
    bind_rows(
      samples_dirichlet_unif |> 
        mutate(distrib = "2_5_3", p = "(1,1,1)/3")
    ) |> 
    bind_rows(
      samples_dirichlet_p |> 
        mutate(distrib = "2_5_3", p = "(3,2,1)/6")
    ) |>
    mutate(
      distrib = factor(
        distrib,
        labels = c(
          "Uniform" = parse(text = latex2exp::TeX("D(1,1,1)")),
          "2_5_3" = parse(text = latex2exp::TeX("$D(2,5,3)$"))
        )
      ),
      p = factor(
        p,
        labels = c(
          "(1,1,1)/3" = parse(text = latex2exp::TeX("$p=(1,1,1)/3$")),
          "(3,2,1)/6" = parse(text = latex2exp::TeX("$p=(3,2,1)/6$"))
        )
      )
    ),
  mapping = aes(x = A, y = B, z = C, color = category)
) +
  geom_point(alpha = .8, size = .5) +
  scale_colour_manual(name = "category", values = colours) +
  facet_grid(p ~ distrib, labeller = label_parsed, switch = "y") +
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

p
  


# Find the intersection point, using a grid----

generate_simplex_grid <- function(resolution = 100) {
  grid <- expand_grid(
    A = seq(0, 1, length.out = resolution),
    B = seq(0, 1, length.out = resolution)
  ) |> 
    mutate(C = 1 - A - B) |> 
    filter(C >= 0)
  
  as_tibble(grid)
}




# KDE for a given category
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
        labels = c(
          "Uniform" = parse(text = latex2exp::TeX("D(1,1,1)")),
          "2_5_3" = parse(text = latex2exp::TeX("$D(2,5,3)$"))
        )
      ),
      p = factor(
        p,
        labels = c(
          "(1,1,1)/3" = parse(text = latex2exp::TeX("$p=(1,1,1)/3$")),
          "(3,2,1)/6" = parse(text = latex2exp::TeX("$p=(3,2,1)/6$"))
        )
      )
    ),
  mapping = aes(x = A, y = B, z = C)
) +
  geom_point(alpha = .8, size = .5, mapping = aes(color = category)) +
  geom_point(
    data = samples_unif_unif$tb_intersection |> 
      mutate(distrib = "Uniform", p = "(1,1,1)/3") |> 
      bind_rows(
        samples_unif_p$tb_intersection |> 
          mutate(distrib = "Uniform", p = "(3,2,1)/6")
      ) |> 
      bind_rows(
        samples_dirichlet_unif$tb_intersection |> 
          mutate(distrib = "2_5_3", p = "(1,1,1)/3")
      ) |> 
      bind_rows(
        samples_dirichlet_p$tb_intersection |> 
          mutate(distrib = "2_5_3", p = "(3,2,1)/6")
      ) |>
      mutate(
        distrib = factor(
          distrib,
          labels = c(
            "Uniform" = parse(text = latex2exp::TeX("D(1,1,1)")),
            "2_5_3" = parse(text = latex2exp::TeX("$D(2,5,3)$"))
          )
        ),
        p = factor(
          p,
          labels = c(
            "(1,1,1)/3" = parse(text = latex2exp::TeX("$p=(1,1,1)/3$")),
            "(3,2,1)/6" = parse(text = latex2exp::TeX("$p=(3,2,1)/6$"))
          )
        )
      )
  ) +
  scale_colour_manual(values = colours) +
  # global_theme() +
  facet_grid(p ~ distrib, labeller = label_parsed, switch = "y") +
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

p

# if (!dir.exists("figs/")) dir.create("figs")
# ggsave(p, file = "figs/baryc-centr-bal.png", width = 5.5, height = 4)

filename <- "baryc-centr-bal"
ggsave(
  p, file = str_c(path, filename, ".pdf"),
  height = 3.3*1.75, width = 3.25*1.75,
  family = font_family,
  device = cairo_pdf
)
# Crop PDF
system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
