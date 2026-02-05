col_group <- c("#00A08A","#F2AD00", "#1b95e0")
col_categ <- c("#56B4E9", "#D55E00", "#CC79A7")
colA <- col_categ[1] ; colB <- col_categ[2] ; colC <- col_categ[3]
colGpe1 <- col_group[2]
colGpe0 <- col_group[1]

library(dplyr)
library(clue)
library(ggtern)
library(compositions)

source("../scripts/utils.R")

# col_categ <- c("#ffdd55","#944edf","#3fb3b2")
col_categ <- c("#56B4E9", "#D55E00", "#CC79A7")
colA <- col_categ[1] ; colB <- col_categ[2] ; colC <- col_categ[3]
# col_group <- c("#1b95e0","darkred")
col_group <- c(colours[["0"]], colours[["1"]])
colGpe1 <- col_group[2]
colGpe0 <- col_group[1]
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

set.seed(1234)
n <- 100
n0 <- n1 <- n
p0 <- c(0.1, 0.5, 0.4)
p1 <- c(0.5, 0.3, 0.2)
# Sample category
x0 <- sample(c(rep("A", p0[1] * n0), rep("B", p0[2] * n0), rep("C", p0[3] * n0)), replace = FALSE)
x1 <- sample(c(rep("A", p1[1] * n1), rep("B", p1[2] * n1), rep("C", p1[3] * n1)), replace = FALSE)
cat_levels <- c("A", "B", "C")

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
tb_sample_z_1

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
tb_sample_z_2

p <- ggtern(
  data = tb_sample_z_1 |> mutate(type = "(1)") |> 
    bind_rows(
      tb_sample_z_2 |> mutate(type = "(2)")
    ), 
  mapping = aes(x = A, y = B, z = C)) +
  geom_point(size = 1, alpha = 0.7, mapping = aes(color = group)) +
  scale_colour_manual(name = "group",values = col_group) +
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

x0_index <- match(x0, cat_levels)
x1_index <- match(x1, cat_levels)
cost_matrix <- outer(x0_index, x1_index, function(i, j) abs(i - j))
# 1-1 matching
assignment <- solve_LSAP(cost_matrix)
# Store this in a tibble
tb_coupling <- tibble(
  x0 = x0,
  x1 = x1[assignment]
) |> 
  mutate(
    cost = abs(match(x0, cat_levels) - match(x1, cat_levels))
  )

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

idx <- which(tb_coupling$cost != 0)

# Ise the plot from previous figure as a baseline
p_matching <- p

# Draw a line joining the matched observations.
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
      color = col_group[1], linewidth = .2,, alpha = .5,
      arrow = arrow(length = unit(0.20, "cm"))
    )
}

p_matching

# filename <- "ternary-categ-ot"
# ggsave(
#   p_matching, file = str_c(path, filename, ".pdf"),
#   height = 2*1.75, width = 3.75*1.75,
#   family = font_family,
#   device = cairo_pdf
# )
# # Crop PDF
# system(paste0("pdfcrop ", path, filename, ".pdf ", path, filename, ".pdf"))
