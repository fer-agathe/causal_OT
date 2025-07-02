library(dplyr)
library(stringr)
library(ggtern)
library(compositions)

# Counterfactuals for categorical features

couleur = c("#1b95e0","darkred")
font_size <- 20
font_family <- "CMU Serif"

path <- "./figs/"
if (!dir.exists(path)) dir.create(path)

theme_ggtern_paper <- function(...) {
  font_family <- "CMU Serif"
  font_size <- 20
  theme(
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text.x = element_text(colour = "black"),
    strip.text = ggtext::element_markdown(),
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    axis.title = element_text(size = rel(.8)),
    tern.axis.arrow.show = TRUE,
    tern.axis.arrow.sep = .13,
    tern.axis.vshift = .05,
    panel.border = element_rect(colour = NA)
  )
}

theme_paper <- function(...) {
  font_family <- "CMU Serif"
  font_size <- 20
  theme(
    text = element_text(family = font_family, size = unit(font_size, "pt")),
    plot.background = element_rect(fill = "transparent", color = NA),
    # legend.text = element_text(size = rel(.8)),
    # legend.title = element_text(size = rel(.8)),
    # legend.title = element_text(size = .8*font_size),
    # legend.text = element_text(size = .8*font_size)
    legend.key = element_blank(),
    panel.spacing = unit(1, "lines"),
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey80"),
    plot.title = element_text(hjust = 0, size = rel(1.3), face = "bold"),
    plot.title.position = "plot",
    strip.background = element_rect(fill = NA, colour = NA)
    # strip.text = element_text(size = rel(1))
  )
}

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
# Figure 1----

## Left Panel----
# Alluvial plot
library(ggalluvial)

set.seed(1234)
n <- 100
n0 <- n1 <- n
p0 <- c(0.5, 0.3, 0.2)
p1 <- c(0.1, 0.5, 0.4)
# Sample category
x0 <- sample(c(rep("A", p0[1] * n0), rep("B", p0[2] * n0), rep("C", p0[3] * n0)), replace = FALSE)
x1 <- sample(c(rep("A", p1[1] * n1), rep("B", p1[2] * n1), rep("C", p1[3] * n1)), replace = FALSE)
cat_levels <- c("A", "B", "C")


# Compute distance between obs from group 0 to group 1
# Setting numeric values to each category: A=1, B=2, C=3
x0_index <- match(x0, cat_levels)
x1_index <- match(x1, cat_levels)
cost_matrix <- outer(x0_index, x1_index, function(i, j) abs(i - j))

min_n <- min(n0, n1)
cost_matrix_square <- cost_matrix[, 1:min_n]

# 1-1 matching
library(clue)
assignment <- solve_LSAP(cost_matrix_square)
# Store this in a tibble
tb_coupling <- tibble(
  x0 = x0,
  x1 = x1[assignment]
) |> 
  mutate(
    cost = abs(match(x0, cat_levels) - match(x1, cat_levels))
  )

tb_coupling |> 
  group_by(x1, x0) |> 
  count() |> 
  group_by(x1) |> 
  mutate(prop = 100 * n / sum(n))


colours_categ <- c("#ffdd55","#944edf","#3fb3b2")

flow <- data.frame(
  depart = rep(LETTERS[1:3], 3),
  category = rep(LETTERS[1:3], each = 3),
  freq = as.vector(table(tb_coupling$x1, tb_coupling$x0))
)

p_1 <- ggplot(
  data = flow,
  mapping = aes(axis1 = depart, axis2 = category, y = freq)
) +
  geom_alluvium(aes(fill = category)) +
  geom_stratum() +
  scale_fill_manual(values = colours_categ) +
  geom_text(
    stat = "stratum",
    mapping = aes(label = after_stat(stratum)),
    family = font_family
  ) +
  # scale_x_discrete(
  #   limits = c("Group 1", "Group 0"),
  #   expand = c(0.15, 0.05)
  # ) +
  scale_x_discrete(
    limits = c("Group 1", "Group 0"),
    labels = c(
      "Group 1" = str_c("<span style='color:", couleur[2], ";'>Group 1</span>"),
      "Group 0" = str_c("<span style='color:", couleur[1], ";'>Group 0</span>")
    ),
    expand = c(0.15, 0.05)
  ) +
  ylab("proportions") + 
  scale_y_continuous(transform = ) +
  # theme_minimal(base_size = font_size, base_family = font_family) +
  theme_paper() +
  theme(
    axis.text.x = ggtext::element_markdown()
  )
p_1


## Right panel----

# dummy dataset to create an empty ternary plot
SB <- tibble(
  A = c(0.2, 0.3, 0.5, 0.6),
  B = c(0.3, 0.4, 0.2, 0.1),
  C = 1 - c(0.2, 0.3, 0.5, 0.6) - c(0.3, 0.4, 0.2, 0.1),
  group = c("1", "1", "0", "0")
)

p_2 <- ggtern(data = SB, aes(x = A, y = B, z = C)) +
  # fake (invisible) points
  # geom_point(size = 0.01, alpha = 0, aes(color = group)) +
  # fake (invisible) lines
  # geom_path(aes(color = group), data = SB, alpha = 0, show.legend = TRUE) +
  scale_colour_manual(name = "group", values = couleur) +
  guides(
    colour = guide_legend(
      override.aes = list(
        linetype = "solid",
        shape = NA,
        size = 1.5,
        alpha = 1
      )
    )
  ) +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_paper() +
  theme(
    legend.title = element_text(size = font_size),
    legend.text = element_text(size = font_size)
    # tern.axis.hshift = .10
  ) +
  theme_latex(TRUE) +
  theme_hidetitles()


p_2 <- p_2 + 
  geom_text(mapping = aes(x = 0.9, y = 0.06, z = 0.08), label = p0[1], color = couleur[1], family = font_family, size = font_size-3, size.unit = "pt") +
  geom_text(mapping = aes(x = 0.09, y = 0.9, z = 0.09), label = p0[2], color = couleur[1], family = font_family, size = font_size-3, size.unit = "pt") +
  geom_text(mapping = aes(x = 0.08, y = 0.06, z = 0.9), label = p0[3], color = couleur[1], family = font_family, size = font_size-3, size.unit = "pt") + 
  geom_text(mapping = aes(x = 0.3, y = 0.1, z = 0.11), label = p1[1], color = couleur[2], family = font_family, size = font_size-3, size.unit = "pt") +
  geom_text(mapping = aes(x = 0.15, y = 0.65, z = 0.25), label = p1[2], color = couleur[2], family = font_family, size = font_size-3, size.unit = "pt") +
  geom_text(mapping = aes(x = 0.1, y = 0.2, z = 0.8), label = p1[3], color = couleur[2], family = font_family, size = font_size-3, size.unit = "pt") 


Li1 <- f_line_simplex(x = c(.75, .125, .125), y = c(.125, .125, .75), lgt = 2)
Li2 <- f_line_simplex(x = c(.75, .125, .125), y = c(.125, .75, .125), lgt = 2)
p_2 <- p_2 + 
  geom_line(
    data = Li2, aes(x = A, y = B, z = C), 
    color = couleur[2], linwidth = .6,
    arrow = arrow(length=unit(0.20,"cm"))
  ) + 
  geom_line(
    data = Li1, aes(x = A, y = B, z = C), 
    color = couleur[2], linwidth = .6,
    arrow = arrow(length=unit(0.20,"cm"))
  ) 
p_2

p_matching_indiv <- cowplot::plot_grid(
  ggplotGrob(
    p_1 +
      # Remove top/bottom margin
      theme(
        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
      )
  ),
  # table_grob,
  ggplotGrob(
    p_2 +
      # Remove top/bottom margin
      theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 0)
      )
  ),
  rel_widths = c(1.4,1),
  ncol = 2
)

p_matching_indiv

filename <- "ternary-categ-matching-indiv"
ggsave(
  p_matching_indiv, file = str_c(path, filename, ".pdf"),
  height = 2*1.75, width = 3.75*1.75,
  family = font_family,
  device = cairo_pdf
)
# Crop PDF
system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))



# Generate data----
set.seed(1234)
n <- 100
n0 <- n1 <- n
p0 <- c(0.5, 0.3, 0.2)
p1 <- c(0.1, 0.5, 0.4)
# Sample category
x0 <- sample(c("A", "B", "C"), size = n0, replace = TRUE, prob = p0)
x1 <- sample(c("A", "B", "C"), size = n1, replace = TRUE, prob = p1)
cat_levels <- c("A", "B", "C")
# Draw data from Dirichlet distributions
library(MCMCpack)

## First type
set.seed(12345)
alpha_A <- c(9, 3, 2)
Z_A <- as.data.frame(rdirichlet(n0 + n1, alpha_A))
alpha_B <- c(3, 11, 4)
Z_B <- as.data.frame(rdirichlet(n0 + n1, alpha_B))
alpha_C <- c(2,3,9)
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


## Second type
set.seed(1234)
alpha_A <- c(19, 3, 2)
Z_A <- as.data.frame(rdirichlet(n0 + n1, alpha_A))
alpha_B <- c(3, 17, 2)
Z_B <- as.data.frame(rdirichlet(n0 + n1, alpha_B))
alpha_C <- c(2, 3, 17)
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

tb_sample_z_2 <- tb_sample_z
rm(tb_sample_z)

p <- ggtern(
  data = tb_sample_z_1 |> mutate(type = "(1)") |> 
    bind_rows(
      tb_sample_z_2 |> mutate(type = "(2)")
    ), 
  mapping = aes(x = A, y = B, z = C)) +
  geom_point(size = 1, alpha = 0.7, mapping = aes(color = group)) +
  scale_colour_manual(name = "group",values = couleur) +
  facet_wrap(~ type) +
  theme_light(base_size = font_size, base_family = font_family) +
  theme_ggtern_paper() +
  theme(
    legend.title = element_text(size = .8 * font_size),
    legend.text = element_text(size = .8 * font_size),
    tern.axis.vshift = .08,
    tern.axis.arrow.sep = .16,
  ) +
  # theme_latex(TRUE)
  theme_hidetitles()
p
filename <- "ternary-categ-drawn"
ggsave(
  p, file = str_c(path, filename, ".pdf"),
  height = 2*1.75, width = 3.75*1.75,
  family = font_family,
  device = cairo_pdf
)
# Crop PDF
system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))

# Matching----

library(compositions)

# Compute distance between obs from group 0 to group 1
# Setting numeric values to each category: A=1, B=2, C=3
x0_index <- match(x0, cat_levels)
x1_index <- match(x1, cat_levels)
cost_matrix <- outer(x0_index, x1_index, function(i, j) abs(i - j))

min_n <- min(n0, n1)
cost_matrix_square <- cost_matrix[, 1:min_n]

# 1-1 matching
library(clue)
assignment <- solve_LSAP(cost_matrix_square)
# Store this in a tibble
tb_coupling <- tibble(
  x0 = x0,
  x1 = x1[assignment]
) |> 
  mutate(
    cost = abs(match(x0, cat_levels) - match(x1, cat_levels))
  )

p_matching <- p

idx <- which(tb_coupling$cost != 0)
for (i in idx) {
  lines_1 <- f_line_simplex(
    x = tb_sample_z_1[i, 1:3], 
    y = tb_sample_z_1[n + assignment[i], 1:3], 
    lgt = 101
  )
  lines_2 <- f_line_simplex(
    x = tb_sample_z_2[i, 1:3], 
    y = tb_sample_z_2[n + assignment[i], 1:3], 
    lgt = 101
  )
  lines_both <- as_tibble(lines_1) |> mutate(type = "(1)") |> 
    bind_rows(
      as_tibble(lines_2) |> mutate(type = "(2)")
    )
  
  p_matching <- p_matching + 
    geom_line(
      data = lines_both, 
      mapping = aes(x = A, y = B, z = C), 
      color = couleur[2], linewidth = .2,, alpha = .5,
      arrow = arrow(length = unit(0.20, "cm"))
    )
}

p_matching

filename <- "ternary-categ-ot"
ggsave(
  p_matching, file = str_c(path, filename, ".pdf"),
  height = 2*1.75, width = 3.75*1.75,
  family = font_family,
  device = cairo_pdf
)
# Crop PDF
system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))


