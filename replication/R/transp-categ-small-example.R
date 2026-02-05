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

theme_ggtern <- function(...) {
  font_family <- "CMU Serif"
  font_size <- 20
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
library(clue)

group_0 <- tribble(
  ~i, ~x, ~p_A, ~p_B, ~p_C, ~y,
  1, "A", .8,   .1,   .1,   1,
  2, "A", .7,   .2,   .1,   2,
  3, "A", .6,   .1,   .3,   3,
  4, "B", .2,   .7,   .1,   4,
  5, "B", .3,   .6,   .1,   5,
  6, "C", .1,   .2,   .7,   6
)

group_1 <- tribble(
  ~i, ~x, ~p_A, ~p_B, ~p_C, ~y,
  7,  "A", .7,  .2,   .1,   3,
  8,  "B", .1,  .7,   .2,   4,
  9,  "B", .6,  .2,   .2,   5,
  10, "B", .2,  .5,   .3,   6,
  11, "C", .3,  .1,   .6,   7,
  12, "C", .1,  .3,   .6,   8
)

matched_ex_1 <- tribble(
  ~i_0, ~i_1,
  1, 12,
  2, 9,
  3, 7,
  4, 10,
  5, 8,
  6, 11
) |>
  left_join(
    group_0 |>
      rename_with(~str_c(.x, "_0")),
    by = "i_0"
  ) |>
  left_join(
    group_1 |>
      rename_with(~str_c(.x, "_1")),
    by = "i_1"
  )

matched_ex_2 <- tribble(
  ~i_0, ~i_1,
  1, 8,
  2, 7,
  3, 11,
  4, 10,
  5, 9,
  6, 12
) |>
  left_join(
    group_0 |>
      rename_with(~str_c(.x, "_0")),
    by = "i_0"
  ) |>
  left_join(
    group_1 |>
      rename_with(~str_c(.x, "_1")),
    by = "i_1"
  )

cat_levels <- c("A", "B", "C")
x0_index <- match(group_0$x, cat_levels)
x1_index <- match(group_1$x, cat_levels)

cost_matrix <- outer(x0_index, x1_index, function(i, j) sqrt((i - j)^2))
cost_matrix

assignment <- solve_LSAP(cost_matrix)

n0 <- nrow(group_0)
n1 <- nrow(group_1)
ot_plan_num <- tibble(
  from = 1:n0,
  to = as.numeric(assignment)
)
ot_plan_num$i_0 <- group_0$i[ot_plan_num$from]
ot_plan_num$i_1 <- group_1$i[ot_plan_num$to]
ot_plan_num

ggtern(
  data = group_0 |> mutate(group = "0") |> 
    bind_rows(group_1 |> mutate(group = "1")) |> 
    left_join(
      ot_plan_num |> 
        mutate(
          id_match = as.character(row_number())
        ) |> 
        dplyr::select(i_0, i_1, id_match) |> 
        pivot_longer(cols = c(i_0, i_1), values_to = "i") |> 
        dplyr::select(-name),
      by = "i"
    ),
  mapping = aes(x = p_A, y = p_B, z = p_C, group = id_match)
) +
  geom_point(mapping = aes(shape = group, colour = x), size = 4) +
  geom_text(
    mapping = aes(
      label = i, 
      x = p_A + ifelse(group == 0, -1, 1) * 0.05
    ),
    size = .3*font_size
  ) +
  labs(x = "$p_A$", y = "$p_B$", z = "$p_C") +
  geom_line(
    colour = "gray40", 
    # mapping = aes(linetype = id_match)
  ) +
  scale_colour_manual(
    name = "category", 
    values = col_categ, 
    labels = c("A" = "A=1", "B" = "B=2", "C" = "C=3")
  ) +
  scale_shape_discrete(name = "group") +
  theme_light(base_size = font_size, base_family = font_family) +
  # theme_paper() +
  theme_ggtern() +
  theme(
    legend.title = element_text(size = .8*font_size),
    legend.text = element_text(size = .8*font_size)
  ) +
  theme_latex(TRUE) +
  theme_hidetitles()

library(compositions)
all_coords <- rbind(
  as.matrix(group_0[, c("p_A", "p_B", "p_C")]),
  as.matrix(group_1[, c("p_A", "p_B", "p_C")])
) |> 
  clr()

row.names(all_coords) <- c(group_0$i, group_1$i)
# Euclidean distances between the clr transform of the propensities
D <- as.matrix(dist(all_coords, method = "euclidean"))
n0 <- nrow(group_0)
n1 <- nrow(group_1)
between_distances <- D[1:n0, (n0 + 1):(n0 + n1)]
round(between_distances, 2)

# source weights
mass_source <- rep(1 / n0, n0)
# target weights
mass_target <- rep(1 / n1, n1)

# Solve the optimal transport plan
ot_plan <- transport::transport(
  a = mass_source, b = mass_target, costm = between_distances, 
  method = "networkflow"
)
ot_plan$i_0 <- group_0$i[ot_plan$from]
ot_plan$i_1 <- group_1$i[ot_plan$to]
ot_plan

all_data <- group_0 |> mutate(group = "0") |> 
  bind_rows(group_1 |> mutate(group = "1"))

library(ggtern)

p <- ggtern(
  data = all_data |> 
    left_join(
      ot_plan |> 
        mutate(
          i_0 = as.numeric(i_0),
          i_1 = as.numeric(i_1),
          id_match = as.character(row_number())
        ) |> 
        dplyr::select(i_0, i_1, id_match) |> 
        pivot_longer(cols = c(i_0, i_1), values_to = "i") |> 
        dplyr::select(-name)
    ),
  mapping = aes(x = p_A, y = p_B, z = p_C, group = id_match)
) +
  geom_point(mapping = aes(shape = group, colour = x), size = 4) +
  geom_text(
    mapping = aes(
      label = i, 
      x = p_A + ifelse(group == 0, -1, 1) * 0.05
    ),
    size = .3*font_size
  ) +
  labs(x = "$p_A$", y = "$p_B$", z = "$p_C") +
  geom_line(
    colour = "gray40", 
    # mapping = aes(linetype = id_match)
  ) +
  scale_colour_manual(name = "category", values = col_categ) +
  scale_shape_discrete(name = "group") +
  theme_light(base_size = font_size, base_family = font_family) +
  # theme_paper() +
  theme_ggtern() +
  theme(
    legend.title = element_text(size = .8*font_size),
    legend.text = element_text(size = .8*font_size)
  ) +
  theme_latex(TRUE) +
  theme_hidetitles()

p

filename <- "ternary-toy"
ggsave(
  p, file = str_c(path, filename, ".pdf"),
  height = 2.2*1.75, width = 4*1.75,
  family = font_family,
  device = cairo_pdf
)
# Crop PDF
system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))

p0 <- as.matrix(group_0[, c("p_A", "p_B", "p_C")])
p1 <- as.matrix(group_1[, c("p_A", "p_B", "p_C")])
p0 <- pmax(p0, 1e-10)
p1 <- pmax(p1, 1e-10)

cross_entropy_cost <- function(x, y) {
  d <- length(x)
  log(mean(y / x)) - mean(log(y / x))
}

between_distances_ce <- outer(
  1:nrow(p0), 1:nrow(p1),
  Vectorize(function(i, j) cross_entropy_cost(p0[i, ], p1[j, ]))
)
round(between_distances_ce, 4)

ot_plan_ce <- transport::transport(
  a = mass_source,
  b = mass_target,
  costm = between_distances_ce,
  method = "networkflow"
)
ot_plan_ce$i_0 <- group_0$i[ot_plan_ce$from]
ot_plan_ce$i_1 <- group_1$i[ot_plan_ce$to]
ot_plan_ce

ggtern(
  data = all_data |> 
    left_join(
      ot_plan_ce |> 
        mutate(
          i_0 = as.numeric(i_0),
          i_1 = as.numeric(i_1),
          id_match = as.character(row_number())
        ) |> 
        dplyr::select(i_0, i_1, id_match) |> 
        pivot_longer(cols = c(i_0, i_1), values_to = "i") |> 
        dplyr::select(-name)
    ),
  mapping = aes(x = p_A, y = p_B, z = p_C, group = id_match)
) +
  geom_point(mapping = aes(shape = group, colour = x), size = 4) +
  geom_text(
    mapping = aes(
      label = i, 
      x = p_A + ifelse(group == 0, -1, 1) * 0.05
    ),
    size = .3*font_size
  ) +
  labs(x = "$p_A$", y = "$p_B$", z = "$p_C") +
  geom_line(
    colour = "gray40", 
    # mapping = aes(linetype = id_match)
  ) +
  scale_colour_manual(name = "category", values = col_categ) +
  scale_shape_discrete(name = "group") +
  theme_light(base_size = font_size, base_family = font_family) +
  # theme_paper() +
  theme_ggtern() +
  theme(
    legend.title = element_text(size = .8*font_size),
    legend.text = element_text(size = .8*font_size)
  ) +
  theme_latex(TRUE) +
  theme_hidetitles()
