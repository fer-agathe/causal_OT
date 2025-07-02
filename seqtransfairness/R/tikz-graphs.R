#' Generates the coordinates used to plot the causal graph from an adjacency
#' matrix
#'
#' @param adj_mat Adjacency matrix, with named columns and rows.
#' @returns A list with three elements:
#' * `coordinates`: Cordinates of the nodes.
#' * `edges`: Edges.
#' * `name_nodes`: Names of the nodes.
#' @importFrom igraph graph_from_adjacency_matrix layout_with_fr
#' @importFrom stringr str_replace_all
generate_coordinates_graph <- function(adj_mat) {
  G <- igraph::graph_from_adjacency_matrix(adj_mat)
  L <- igraph::layout_with_fr(G)
  N <- rownames(adj_mat)
  N <- str_replace_all(N, "_", " ")
  E <- matrix(NA, sum(adj_mat), 2)
  i_s = 0
  for (i in 1:nrow(adj_mat)) {
    for (j in which(adj_mat[i, ] > 0)) {
      i_s <- i_s + 1
      E[i_s, ] <- c(i, j)
    }
  }
  edges <- list()
  for(i in 1:nrow(E)) edges[[i]] <- E[i,]

  CRDNT <- L
  CRDNT[, 1] <- (L[, 1] - min(L[, 1])) * 10 / (max(L[, 1]) - min(L[, 1]))
  CRDNT[, 2] <- (L[, 2] - min(L[, 2])) * 5 / (max(L[, 2]) - min(L[, 2]))

  coordinates <- list()
  for (i in 1:nrow(CRDNT)) coordinates[[i]] <- CRDNT[i, ]
  name_nodes <- N
  list(
    coordinates = coordinates,
    edges = edges,
    name_nodes = name_nodes
  )
}

#' Function to generate TikZ code for nodes
#'
#' @param coordinates Coordinates of the nodes.
#' @param name_nodes Vector of names of the nodes.
#' @param colour_nodes Colour of nodes. Default to "yellow!60".
generate_nodes <- function(coordinates,
                           name_nodes,
                           colour_nodes = "yellow!60") {
  tikz_nodes <- ""
  if (length(colour_nodes) == 1)
    colour_nodes <- rep(colour_nodes, length(name_nodes))
  for (i in seq_along(coordinates)) {
    x <- coordinates[[i]][1]
    y <- coordinates[[i]][2]
    nm <- name_nodes[i]
    cn <- colour_nodes[i]
    ligne <- paste0(
      "\\node[fill=", cn,
      "]",
      sprintf(" (n%d) at (%.1f, %.1f) {", i, x, y),
      nm, "};\n"
    )
    tikz_nodes <- paste0(tikz_nodes, ligne)
  }
  tikz_nodes
}

#' Function to generate TikZ code for edges
#'
#' @param edges Index of edges for each node.
generate_edges <- function(edges) {
  tikz_edges <- ""
  for (edge in edges) {
    tikz_edges <- paste0(
      tikz_edges,
      sprintf("\\draw[->, black] (n%d) -- (n%d);\n", edge[1], edge[2])
    )
  }
  tikz_edges
}

#' Tikz for Causal Graph
#'
#' @param adj_mat Adjacency matrix, with named columns and rows.
#' @param colour_nodes Colour of nodes. Default to "yellow!60". If specific
#'        colour is desired for a node, a vector with colours in the same order
#'        as the column names must be provided.
#' @export
causal_graph_tikz <- function(adj_mat, colour_nodes = "yellow!60") {
  coords_g <- generate_coordinates_graph(adj_mat)
  coordinates <- coords_g$coordinates
  edges <- coords_g$edges
  name_nodes <- coords_g$name_nodes
  nodes <- generate_nodes(coordinates, name_nodes, colour_nodes)
  edges <- generate_edges(edges)
  tikz_code <- paste0(
    "\\begin{tikzpicture}\n", nodes, edges, "\\end{tikzpicture}"
  )
  tikz_code
}


