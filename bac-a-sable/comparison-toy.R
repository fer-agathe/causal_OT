library(tidyverse)

library(extrafont)
loadfonts(device = "pdf")
colours <- c("A" = "#FFE788", "B" = "#9F78C4", "C" = "#5FC7C6")
font_size <- 20
font_family <- "CMU Serif"

path <- "./figs/"
if (!dir.exists(path)) dir.create(path)

theme_paper <- function(...) {
  font_family <- "CMU Serif"
  font_size <- 20
  theme(
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.text = element_text(size = rel(1.1)),
    legend.title = element_text(size = rel(1.1)),
    #legend.background = element_rect(
    #  fill = "transparent", linetype = "solid", colour = "black"
    #),
    # legend.position = "bottom",
    # legend.direction = "horizontal",
    # legend.box = "vertical",
    legend.key = element_blank(),
    panel.spacing = unit(1, "lines"),
    plot.title = element_text(hjust = 0, size = rel(1.3), face = "bold"),
    plot.title.position = "plot",
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text = element_text(size = rel(1.1))
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


group_0 <- tribble(
  ~i, ~x, ~p_A, ~p_B, ~p_C, ~y,
  1, "A", .8, .1, .1, 1,
  2, "A", .7, .2, .1, 2,
  3, "A", .6, .1, .3, 3,
  4, "B", .2, .7, .1, 4,
  5, "B", .3, .6, .1, 5,
  6, "C", .1, .2, .7, 6
)

group_1 <- tribble(
  ~i, ~x, ~p_A, ~p_B, ~p_C, ~y,
  7,  "A", .7, .2, .1, 3,
  8,  "B", .1, .7, .2, 4,
  9,  "B", .6, .2, .2, 5,
  10, "B", .2, .5, .3, 6,
  11, "C", .3, .1, .6, 7,
  12, "C", .1, .3, .6, 8
)

# First random matching----
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

matched_ex_1 |>
  mutate(diff = y_1 - y_0) |>
  select(i_0, i_1, diff)

# Second random matching----
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

matched_ex_2 |>
  mutate(diff = y_1 - y_0) |>
  select(i_0, i_1, diff)

# Distances----

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
between_distances <- D[(n0 + 1):(n0 + n1), 1:n0]
round(between_distances, 2)


# Matching----

# source weights
mass_source <- rep(1 / n1, n1)
# target weights
mass_target <- rep(1 / n0, n0)

# Solve the optimal transport plan
ot_plan <- transport::transport(
  a = mass_source, b = mass_target, costm = between_distances, 
  method = "networkflow"
)

ot_plan$i_0 <- group_0$i[ot_plan$to]
ot_plan$i_1 <- group_1$i[ot_plan$from]
ot_plan


all_data <- 
  bind_rows(group_0 |> mutate(group = "0"), group_1 |> mutate(group = "1"))

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
  scale_colour_manual(name = "category", values = colours) +
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

ot_plan |> 
  left_join(group_0 |> select(i_0 = i, y_0 = y), by = "i_0") |> 
  left_join(group_1 |> select(i_1 = i, y_1 = y), by = "i_1") |> 
  mutate(diff = y_1 - y_0) |>
  select(i_0, i_1, diff)


