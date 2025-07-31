#' seqtransfairness
#'
#' @name seqtransfairness
#' @description
#' {seqtransfairness} is an R package that provides tools to create
#' counterfactuals based on sequential conditional transport using a causal
#' graph assumed to be known. This approach links two existing approaches
#' to derive counterfactuals: adaptations based on a causal graph,
#' and optimal transport.
#'
#' @author Fernandes Machado, Agathe and Charpentier, Arthur and Gallic, Ewen.
#'
#' @examples
#' # Simulate data
#' sim_dat <- simul_dataset()
#' # Causal graph:
#' variables <- c("S", "X1", "X2", "Y")
#' adj <- matrix(
#'   # S  X1 X2 Y
#'   c(0, 1, 1, 1,# S
#'     0, 0, 1, 1,# X1
#'     0, 0, 0, 1,# X2
#'     0, 0, 0, 0  # Y
#'   ),
#'   ncol = length(variables),
#'   dimnames = rep(list(variables), 2),
#'   byrow = TRUE
#' )
#' # To visualize the causal graph:
#' # causal_graph <- fairadapt::graphModel(adj)
#' # plot(causal_graph)
#'
#' # Sequential transport according to the causal graph
#' transported <- seq_trans(data = sim_dat, adj = adj, s = "S", S_0 = 0, y = "Y")
#' transported
#' # Transported values from S=0 to S=1, using the causal graph.
#' predict(transported)
#' # Transport of two new observations
#' new_obs <- data.frame(X1 = c(-2, -1), X2 = c(-1, -2), S = c(0, 0))
#' predict(transported, newdata = new_obs, data = sim_dat)
NULL

#' Prints the Transported Values of a `sequential_transport` Object
#'
#' @method print sequential_transport
#' @export
#' @noRd
print.sequential_transport <- function(x, ...) {
  params <- x$params
  top_order <- params$top_order
  cat("Outcome:", params$y, "\n")
  cat("Sensitive:", params$s, "\n")
  cat("\t- with value", params$S_0, "in source group.\n")
  cat("\t- twith value", params$S_1, "in target group.\n")
  cat("\nTopological ordering:", top_order, "\n")

  # From fairadapt code
  mat <- rbind(colnames(params$adj), params$adj)
  mat <- apply(mat, 2L, format, justify = "right")
  mat <- cbind(format(c("", rownames(params$adj))), mat)

  cat("\nBased on causal graph:\n")
  cat(apply(mat, 1L, paste, collapse = " "), sep = "\n")
  cat("\n")
  invisible(x)
}


#' Transported Values by Sequential Transport
#'
#' @param object Object of class `"sequential_transport"`.
#' @param newdata An optional data frame/tibble in which to look for variables
#'        with which to predict.
#' @param data Initial data on which the object has been created.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
#' @md
#' @importFrom tibble as_tibble tibble
#' @importFrom dplyr mutate
#' @importFrom rlang !! sym :=
#' @examples
#' sim_dat <- simul_dataset()
#' # Causal graph:
#' variables <- c("S", "X1", "X2", "Y")
#' adj <- matrix(
#'   # S  X1 X2 Y
#'   c(0, 1, 1, 1,# S
#'     0, 0, 1, 1,# X1
#'     0, 0, 0, 1,# X2
#'     0, 0, 0, 0  # Y
#'   ),
#'   ncol = length(variables),
#'   dimnames = rep(list(variables), 2),
#'   byrow = TRUE
#' )
#' # To visualize the causal graph:
#' # causal_graph <- fairadapt::graphModel(adj)
#' # plot(causal_graph)
#'
#' # Sequential transport according to the causal graph
#' transported <- seq_trans(data = sim_dat, adj = adj, s = "S", S_0 = 0, y = "Y")
#' # The transported values for group in S=0:
#' predict(transported)
#' # Transport of new observations from S=0:
#' new_obs <- data.frame(S=0, X1 = c(-1, -1.5), X2 = c(0, -.5))
#' predict(transported, newdata = new_obs, data = sim_dat)
predict.sequential_transport <- function(object,
                                         newdata,
                                         data,
                                         ...) {
  no_data <- (missing(newdata) || is.null(newdata))
  if (no_data) {
    pred <- object$transported
  } else {
    if (missing(data)) {
      stop("Please also provide initial data to the `data` argument")
    }
    pred <- seq_trans_new(x = object, newdata = newdata, data = data)
  }
  as_tibble(pred) |>
    mutate(!!sym(object$params$s) := object$params$S_1)
}
