#' Compute counterfactual predicted probabilities for a single observation
#'
#' @param i Index of the observation from group °
#' @param pred_probs_0 Matrix of predicted probabilities (group 0).
#' @param pred_probs_1 Matrix of predicted probabilities (group 1).
#' @param weights_0 Matrix of intra-group distances for group 0.
#' @param weights_1 Matrix of inter-group distances from group 0 to group 1.
#' @param num_neighbors_q Number of neighbors to use.
#'
#' @return A vector of counterfactual predicted probabilities.
compute_counterfactual_probs <- function(i,
                                         pred_probs_0,
                                         pred_probs_1,
                                         weights_0,
                                         weights_1,
                                         num_neighbors_q) {
  # Identify closest neighbours within the same group
  dist_neigh_0 <- weights_0[i, , drop = FALSE]
  # and among the other group
  dist_neigh_1 <- weights_1[i, , drop = FALSE]
  
  ranks_weights_0 <- order(dist_neigh_0, decreasing = TRUE)[1:num_neighbors_q]
  ranks_weights_1 <- order(dist_neigh_1, decreasing = TRUE)[1:num_neighbors_q]
  i_rank <- which(ranks_weights_0 == i)
  
  W_i <- wasserstein_simplex(
    X = pred_probs_0[ranks_weights_0, ],
    Y = pred_probs_1[ranks_weights_1, ],
    wx = dist_neigh_0[ranks_weights_0],
    wy = dist_neigh_1[ranks_weights_1]
  )
  
  # Most likely match under the transport plan
  pred_probs_1[ranks_weights_1, ][which.max(W_i$plan[i_rank, ]), ]
}


#' @param pred_probs_0 Matrix with predicted probabilities to belong the each 
#'  class of the categorical variable, in group 0 (source).
#' @param pred_probs_1 Matrix with predicted probabilities to belong the each 
#'  class of the categorical variable, in group 1 (target).
#' @param weights_0 Weights corresponding to the distance between observations 
#'  within source group.
#' @param weights_1 Weights corresponding to the distance between observations 
#'  within target group.
#' @param num_neighbors_q Number of neigbors to use for categorical variables.
#'  Default to the min between 50 and the number of observations in the data.
#' @param cl A cluster object, created by package parallel. If `NULL` (default), 
#'   no parallel computing is used to transport categorical data.
#'  
#' @importFrom transportsimplex wasserstein_simplex
ot_simplex_probs <- function(pred_probs_0,
                             pred_probs_1,
                             weights_0,
                             weights_1,
                             num_neighbors_q = NULL,
                             cl = NULL) {
  
  if (is.null(num_neighbors_q)) {
    num_neighbors_q <- min(nrow(pred_probs_0), nrow(pred_probs_1), 50)
  } else {
    num_neighbors_q <- min(nrow(pred_probs_0), nrow(pred_probs_1), num_neighbors_q)
  }
  
  mat_counter_categ <- matrix(
    NA, ncol = ncol(pred_probs_0), nrow = nrow(pred_probs_0)
  )
  
  indices <- seq_len(nrow(pred_probs_0))
  
  if (!is.null(cl)) {
    parallel::clusterExport(
      cl, varlist = c(
        "pred_probs_0", "pred_probs_1", 
        "weights_0", "weights_1", 
        "num_neighbors_q", "compute_counterfactual_probs"
      ),
      envir = environment()
    )
    res <- pbapply::pblapply(
      indices, 
      compute_counterfactual_probs,
      pred_probs_0 = pred_probs_0,
      pred_probs_1 = pred_probs_1,
      weights_0 = weights_0,
      weights_1 = weights_1,
      num_neighbors_q = num_neighbors_q,
      cl = cl
    )
  } else {
    res <- pbapply::pblapply(
      indices, compute_counterfactual_probs,
      pred_probs_0 = pred_probs_0,
      pred_probs_1 = pred_probs_1,
      weights_0 = weights_0,
      weights_1 = weights_1,
      num_neighbors_q = num_neighbors_q
    )
    
  }
  
  do.call(rbind, res)
}

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
  factor(c(1, 2, 4), levels = 1:4, labels = c("A", "B", "C", "D"))
  
  factor(assignments, levels = 1:length(labels), labels = labels)
}

#' Sequential Transport Using a Pre-Defined Causal Graph
#'
#' The sensitive attribute, S, is assumed to be a binary variable with value
#' $S_0$ in the source distribution and $S_1$ in the target distribution.
#'
#' @param data Data frame with the observations.
#' @param adj Adjacency matrix for the causal graph.
#' @param s Name of the sensitive attribute column in the data.
#' @param S_0 Label of the sensitive attribute in the source distribution.
#' @param y Name of the outcome variable in the data.
#' @param num_neighbors Number of neighbors to use in the weighted quantile
#'        estimation. Default to 5.
#' @param num_neighbors_q Number of neigbors to use for categorical variables.
#'  Default to the min between 50 and the number of observations in the data.
#' @param silent If `TRUE`, the messages showing progress in the estimation are
#'        not shown. Default to `silent=FALSE`.
#' @param cl A cluster object, created by package parallel. If `NULL` (default), 
#'   no parallel computing is used to transport categorical data. Otherwise, 
#'   only used to transport categorical data.
#'
#' @returns An element of class `"sequential_transport"` (a list):
#' * `transported`: A named list with the transported values. The names are those of the variables.
#' * `weights`: A list with the weights of each observation in the two groups.
#' * `ecdf`: A list with empirical distribution functions for numerical variables.
#' * `ecdf_values`: A list with the values of the ecdf evaluated for each observation in the source distribution.
#' * `fit_for_categ`: A list with the estimated multinomial models to predict categories using parents characteristics
#' * `params`: A list with some parameters used to transport observations:
#'     * `adj`: Adjacency matrix.
#'     * `top_order`: Topological ordering.
#'     * `s`: Name of the sensitive attribute.
#'     * `S_0`: Label of the sensitive attribute in the source distribution.
#'     * `S_1`: Label of the sensitive attribute in the target distribution.
#'     * `y`: Name of the outcome variable in the data.
#'     * `num_neighbors`: Number of neighbors used when computing quantiles.
#' @md
#' @export
#'
#' @examples
#' # Data with two groups: S=0, S=1, an outcome Y and two covariates X1 and X2
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
#' transported_val <- as.data.frame(transported$transported)
#' head(transported_val)
#' @importFrom stats predict ecdf quantile
#' @importFrom dplyr across filter mutate pull select
#' @importFrom tidyselect where
#' @importFrom rlang sym !! := is_character
#' @importFrom cluster daisy
#' @importFrom Hmisc wtd.quantile
#' @importFrom nnet multinom
#' @importFrom purrr map_chr
#' @seealso [seq_trans_new()], [simul_dataset()]
seq_trans <- function(data,
                      adj,
                      s,
                      S_0,
                      y,
                      num_neighbors = 5,
                      num_neighbors_q = NULL,
                      silent = FALSE,
                      cl = NULL) {
  # Make sure character variables are encoded as factors
  data <-
    data |>
    mutate(across(where(is_character), ~as.factor(.x)))
  
  s_unique <- unique(data[[s]])
  S_1 <- s_unique[s_unique != S_0]
  
  # Topological ordering
  top_order <- seqtransfairness::topological_ordering(adj)
  variables <- top_order[!top_order %in% c(s, y)]
  # Observations in group S_0
  data_0 <- data |> filter(!!sym(s) == !!S_0)
  data_1 <- data |> filter(!!sym(s) != !!S_0)
  
  # Lists where results will be stored
  list_transported <- list()  # Transported values
  list_transported_prob <- list()  # Transported prob. for categ. variables
  list_weights <- list()      # Weights
  list_ecdf <- list()         # Empirical dist. function
  list_ecdf_values <- list()  # Evaluated values of the ecdf
  fit_for_categ <- list()     # Fitted multinomial models for categ. variables
  gower_matrix_all <- NULL    # Distance between observations
  
  for (x_name in variables) {
    if (silent == FALSE) cat("Transporting ", x_name, "\n")
    # Names of the parent variables
    parents <- colnames(adj)[adj[, x_name] == 1]
    # values of current x in each group
    x_S0 <- data_0 |> pull(!!x_name)
    x_S1 <- data_1 |> pull(!!x_name)
    # Check whether X is numeric
    is_x_num <- is.numeric(x_S0)
    # Characteristics of the parent variables (if any)
    parents_characteristics <- data_0 |> select(!!parents, -!!s)
    
    if (length(parents_characteristics) > 0) {
      
      data_0_parents <- data_0 |> select(!!parents) |> select(-!!s)
      data_1_parents <- data_1 |> select(!!parents) |> select(-!!s)
      # Weights in S_0
      weights_S0 <- as.matrix(daisy(data_0_parents, metric = "gower"))
      tot_weights_S0 <- apply(weights_S0, MARGIN = 1, sum)
      # Weights in S_1
      # First, we need to get the transported values for the parents, if necessary
      data_0_parents_t <- data_0_parents #init
      for (parent in parents) {
        # does the parent depend on the sensitive variable
        if (parent %in% names(list_transported)) {
          data_0_parents_t <-
            data_0_parents_t |>
            mutate(!!sym(parent) := list_transported[[parent]])
        }
      }
      # Unfortunately, we will compute a lot of distances not needed
      combined <- rbind(data_0_parents_t, data_1_parents)
      gower_dist <- daisy(combined, metric = "gower")
      gower_matrix <- as.matrix(gower_dist)
      n_0 <- nrow(data_0_parents_t)
      n_1 <- nrow(data_1_parents)
      weights_S1 <- gower_matrix[1:n_0, (n_0 + 1):(n_0 + n_1), drop = FALSE]
      weights_S1 <- weights_S1 + 1e-8
      weights_S1 <- 1 / (weights_S1)^2
      tot_weights_S1 <- apply(weights_S1, MARGIN = 1, sum)
      
      if (is_x_num == TRUE) {
        # Numerical variable to transport
        
        # Empirical distribution function
        f <- rep(NA, length(x_S0))
        for (i in 1:length(x_S0)) {
          f[i] <- weights_S0[i, ] %*% (x_S0 <= x_S0[i]) / tot_weights_S0[i]
        }
        list_ecdf_values[[x_name]] <- f
        f[f==1] <- 1-(1e-8)
        
        # Transported values
        transported <- rep(NA, length(x_S0))
        for (i in 1:length(x_S0)) {
          wts <- weights_S1[i, ]
          wts[-order(wts, decreasing = TRUE)[1:num_neighbors]] <- 0
          transported[i] <- Hmisc::wtd.quantile(
            x = x_S1, weights = weights_S1[i, ], probs = f[i]
          ) |> suppressWarnings()
        }
      } else {
        # X is non numeric and has parents
        x_labels <- data |> pull(!!x_name) |> levels()
        # Estimation of propensity in source group
        fit_categ_0 <- nnet::multinom(
          paste(x_name, "~ ."),
          data = data_0 |> select(-!!y),
          trace = FALSE
        )
        # Estimation of propensity in target group
        fit_categ_1 <- nnet::multinom(
          paste(x_name, "~ ."),
          data = data_1 |> select(-!!y),
          trace = FALSE
        )
        
        # Predictions with these models:
        pred_probs_0 <- predict(fit_categ_0, type = "probs")
        pred_probs_1 <- predict(fit_categ_1, type = "probs")
        
        if (length(x_labels) == 2) {
          # Binary
          # Empirical distribution function
          f <- rep(NA, length(pred_probs_0))
          for (i in 1:length(pred_probs_0)) {
            f[i] <- weights_S0[i, ] %*% (pred_probs_0 <= pred_probs_0[i]) / tot_weights_S0[i]
          }
          list_ecdf_values[[x_name]] <- f
          f[f==1] <- 1-(1e-8)
          
          # Transported values
          pred_probs_0_t <- rep(NA, length(pred_probs_0))
          for (i in 1:length(pred_probs_0)) {
            wts <- weights_S1[i, ]
            wts[-order(wts, decreasing = TRUE)[1:num_neighbors]] <- 0
            pred_probs_0_t[i] <- Hmisc::wtd.quantile(
              x = pred_probs_1, weights = weights_S1[i, ], probs = f[i]
            ) |> suppressWarnings()
          }
          pred_probs_0_t <- cbind(pred_probs_0_t, 1-pred_probs_0_t)
          transported <- get_assignment(
            probs = pred_probs_0_t, 
            labels = x_labels, 
            p = table(data_1 |> pull(!!x_name)) / nrow(data_1)
          )
        } else {
          # Categorical with more than two classes
          
          # If some classes are in a group but not in the other
          if (!all(x_labels %in% colnames(pred_probs_0))) {
            small_prob <- min(pred_probs_0/2, 1e-8)
            # Identify missing columns
            x_labels_missing_0 <- 
              x_labels[which(! x_labels %in% colnames(pred_probs_0))]
            # set those to a tiny value
            pred_probs_0_missing <- matrix(
              rep(small_prob, n_0 * length(x_labels_missing_0)), 
              ncol = length(x_labels_missing_0)
            )
            colnames(pred_probs_0_missing) <- x_labels_missing_0
            # Add column(s) to the prediction matrix
            pred_probs_0 <- cbind(pred_probs_0, pred_probs_0_missing)
            # Normalize the probabilities
            pred_probs_0 <- pred_probs_0 / rowSums(pred_probs_0)
          }
          
          # Same for other group
          if (!all(x_labels %in% colnames(pred_probs_1))) {
            small_prob <- min(pred_probs_1/2, 1e-8)
            # Identify missing columns
            x_labels_missing_1 <- 
              x_labels[which(! x_labels %in% colnames(pred_probs_1))]
            # set those to a tiny value
            pred_probs_1_missing <- matrix(
              rep(small_prob, n_1 * length(x_labels_missing_1)), 
              ncol = length(x_labels_missing_1)
            )
            colnames(pred_probs_1_missing) <- x_labels_missing_1
            # Add column(s) to the prediction matrix
            pred_probs_1 <- cbind(pred_probs_1, pred_probs_1_missing)
            # Normalize the probabilities
            pred_probs_1 <- pred_probs_1 / rowSums(pred_probs_1)
          }
          
          pred_probs_0_t <- ot_simplex_probs(
            pred_probs_0 = pred_probs_0, 
            pred_probs_1 = pred_probs_1, 
            weights_0 = 1 / (weights_S0 + 1e-8)^2, 
            weights_1 = weights_S1, 
            num_neighbors_q = num_neighbors_q,
            cl = cl
          )
          
          # Target prob
          target_prob <- table(data_1 |> pull(!!x_name)) / nrow(data_1)
          
          transported <- get_assignment(
            probs = pred_probs_0_t, 
            labels = x_labels, 
            p = target_prob
          )
        }
        
        fit_for_categ[[x_name]] <- list(
          "source" = fit_categ_0,
          "target" = fit_categ_1
        )
        colnames(pred_probs_0_t) <- x_labels
        list_transported_prob[[x_name]] <- pred_probs_0_t
      }
      list_transported[[x_name]] <- transported
      
      # Store weights for possible later use
      list_weights[[x_name]] <- list(
        w_S0 = list(weights = weights_S0, tot_weights = tot_weights_S0),
        w_S1 = list(weights = weights_S1, tot_weights = tot_weights_S1)
      )
    } else {
      # No parents
      if (is_x_num == TRUE) {
        # X is numerical and has no parents
        F_X_S0 <- ecdf(x_S0)
        list_ecdf[[x_name]] <- F_X_S0
        f <- F_X_S0(x_S0)
        list_ecdf_values[[x_name]] <- f
        transported <- as.numeric(quantile(x_S1, probs = f))
      } else {
        # X is not numerical and has no parents
        x_labels <- data |> pull(!!x_name) |> levels()
        # Estimation of propensity in source group
        fit_categ_0 <- nnet::multinom(
          paste(x_name, "~ ."),
          data = data_0 |> select(-!!y),
          trace = FALSE
        )
        # Estimation of propensity in target group
        fit_categ_1 <- nnet::multinom(
          paste(x_name, "~ ."),
          data = data_1 |> select(-!!y),
          trace = FALSE
        )
        # Predictions with these models:
        pred_probs_0 <- predict(fit_categ_0, type = "probs")
        pred_probs_1 <- predict(fit_categ_1, type = "probs")
        
        if (length(x_labels) == 2) {
          # Binary variable
          F_pred_probs_S0 <- ecdf(pred_probs_0)
          list_ecdf[[x_name]] <- F_pred_probs_S0
          f <- F_pred_probs_S0(pred_probs_0)
          list_ecdf_values[[x_name]] <- f
          pred_probs_0_t <- as.numeric(quantile(pred_probs_1, probs = f))
          pred_probs_0_t <- cbind(pred_probs_0_t, 1-pred_probs_0_t)
          transported <- get_assignment(
            probs = pred_probs_0_t, 
            labels = x_labels, 
            p = table(data_1 |> pull(!!x_name)) / nrow(data_1)
          )
          colnames(pred_probs_0_t) <- x_labels
        } else {
          # More than two classes
          
          # If some classes are in a group but not in the other
          if (!all(x_labels %in% colnames(pred_probs_0))) {
            small_prob <- min(pred_probs_0/2, 1e-8)
            # Identify missing columns
            x_labels_missing_0 <- 
              x_labels[which(! x_labels %in% colnames(pred_probs_0))]
            # set those to a tiny value
            pred_probs_0_missing <- matrix(
              rep(small_prob, n_0 * length(x_labels_missing_0)), 
              ncol = length(x_labels_missing_0)
            )
            colnames(pred_probs_0_missing) <- x_labels_missing_0
            # Add column(s) to the prediction matrix
            pred_probs_0 <- cbind(pred_probs_0, pred_probs_0_missing)
            # Normalize the probabilities
            pred_probs_0 <- pred_probs_0 / rowSums(pred_probs_0)
          }
          
          # Same for other group
          if (!all(x_labels %in% colnames(pred_probs_1))) {
            small_prob <- min(pred_probs_1/2, 1e-8)
            # Identify missing columns
            x_labels_missing_1 <- 
              x_labels[which(! x_labels %in% colnames(pred_probs_1))]
            # set those to a tiny value
            pred_probs_1_missing <- matrix(
              rep(small_prob, n_1 * length(x_labels_missing_1)), 
              ncol = length(x_labels_missing_1)
            )
            colnames(pred_probs_1_missing) <- x_labels_missing_1
            # Add column(s) to the prediction matrix
            pred_probs_1 <- cbind(pred_probs_1, pred_probs_1_missing)
            # Normalize the probabilities
            pred_probs_1 <- pred_probs_1 / rowSums(pred_probs_1)
          }
          
          mapping <- wasserstein_simplex(
            X = pred_probs_0, Y = pred_probs_1
          )
          pred_probs_0_t <- counterfactual_w(
            mapping = mapping, X0 = pred_probs_0, X1 = pred_probs_1
          )
          transported <- get_assignment(
            probs = pred_probs_0_t, 
            labels = x_labels, 
            p = table(data_1 |> pull(!!x_name)) / nrow(data_1)
          )
          colnames(pred_probs_0_t) <- x_labels
        }
        list_transported_prob[[x_name]] <- pred_probs_0_t
      }
      list_transported[[x_name]] <- transported
    }
  }
  
  structure(
    list(
      transported = list_transported,
      transported_prob = list_transported_prob,
      weights = list_weights,
      ecdf = list_ecdf,
      ecdf_values = list_ecdf_values,
      fit_for_categ = fit_for_categ,
      params = list(
        adj = adj,
        top_order = top_order,
        s = s,
        S_0 = S_0,
        S_1 = S_1,
        y = y,
        num_neighbors = num_neighbors
      )
    ),
    class = "sequential_transport"
  )
}

#' Pairwise distance matrix on the simplex
#'
#' @description
#' Computes the pairwise distance matrix of observations in the simplex, using
#' the cost function for optimal transport on the unit simplex as the distance
#' metric.
#'
#' @param X Matrice of observations (one observation per row).
#' @param Y Matrice of observations (one observation per row).
#'
#' @returns A matrix of size m x n, where m is the number of observation in X,
#'  and n is the number of observations in X, containing the distances between
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
#' Finds the optimal transport plan using linear programming.
#' In a first attempts, it uses `CVXR::solve` with the OSQP solver.
#' If this fails, it uses `lpSolve::lp` instead.
#' The function minimizes the transport cost while ensuring:
#' * Mass conservation (row and column sums match the marginals).
#' * Nonnegative transport flows.
#'
#' @param dxy Cost matrix of transport distances between points in X and Y.
#' @param wx Weights (marginal distribution) for X.
#' @param wy Weights (marginal distribution) for Y.
#' @param p Order of the Wassterstein distance. (If p=2: squared Euclidean
#'  cost).
#'
#' @importFrom CVXR Variable Minimize matrix_trace Problem solve
#' @importFrom lpSolve lp
#'
#' @noRd
wass_lp <- function(dxy,
                    wx,
                    wy,
                    p) {
  cxy    <- dxy
  m      <- length(wx)
  ww_m   <- matrix(wx, ncol = 1)
  n      <- length(wy)
  ww_n   <- matrix(wy, nrow = 1)
  ones_m <- matrix(rep(1, n), ncol = 1)
  ones_n <- matrix(rep(1, m), nrow = 1)
  plan   <- CVXR::Variable(m, n)
  
  wd.obj    <- CVXR::Minimize(CVXR::matrix_trace(t(cxy) %*% plan))
  wd.const1 <- list(plan >= 0)
  wd.const2 <- list(plan %*% ones_m == ww_m, ones_n %*% plan == ww_n)
  wd.prob   <- CVXR::Problem(wd.obj, c(wd.const1, wd.const2))
  wd.solve  <- CVXR::solve(wd.prob, solver = "OSQP")
  
  if (all(wd.solve$status=="optimal")) {
    # successful
    gamma <- wd.solve$getValue(plan)
    value <- (base::sum(gamma * cxy))
  } else {
    # failed : use lpsolve
    cxy <- (dxy)
    m   <- nrow(cxy)
    n   <- ncol(cxy)
    
    c  <- as.vector(cxy)
    A1 <- base::kronecker(matrix(1, nrow = 1, ncol = n), diag(m))
    A2 <- base::kronecker(diag(n), matrix(1, nrow = 1, ncol = m))
    A  <- rbind(A1, A2)
    
    f.obj <- c
    f.con <- A
    f.dir <- rep("==", nrow(A))
    f.rhs <- c(rep(1 / m, m), rep(1 / n, n))
    f.sol <- (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
    
    gamma <- matrix(f.sol$solution, nrow = m)
    value <- (sum(gamma*cxy)^(1 / p))
  }
  list(distance = value, plan = gamma)
}


#' Ensures that a weight vector (marginal distribution) is valid
#'
#' @description
#' Returns a uniform weight if the provided vector if NULL. Otherwise, checks
#' if the vector has length M and nonnegative entries, and if so, normalizes
#' the vector of weights to sum to 1.
#'
#' @param mvec (Optional) Vector of weights.
#' @param M Length of the weight vector.
#' @param fname Name of the distance used (string).
#' @noRd
valid_single_marginal <- function(mvec, M, fname) {
  dname <- paste0("'", deparse(substitute(mvec)), "'")
  if ((length(mvec) == 0) && is.null(mvec)) {
    return(rep(1 / M, M))
  } else {
    mvec <- as.vector(mvec)
    if ((length(mvec) != M) || (any(mvec < 0))) {
      stop(
        paste0(
          "* ", fname, " : ", dname,
          " should be a nonnegative vector of length ",M,"."
        )
      )
    }
    return(mvec / base::sum(mvec))
  }
}


#' Wasserstein distance between two sets of probability vectors X and Y
#'
#' @param X Matrix of probability vectors in a first group.
#' @param Y Matrix of probability vectors in a second group.
#' @param wx Weights (marginal distribution) for X. Default to `NULL` (uniform
#' weights will be used).
#' @param wy Weights (marginal distribution) for Y. Default to `NULL` (uniform
#' weights will be used).
#'
#' @returns A list with two elements:
#' * `distance`: the Wassterstein distance
#' * `plan`: the optimal transport plan describing how mass is transported
#'   between X and Y.
#' @export
wasserstein_simplex <- function(X,
                                Y,
                                wx = NULL,
                                wy = NULL) {
  ## CHECK INPUTS
  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }
  if (is.vector(Y)) {
    Y <- matrix(Y, ncol = 1)
  }
  if (!is.matrix(X)) { stop("* wasserstein : input 'X' should be a matrix.") }
  if (!is.matrix(Y)) { stop("* wasserstein : input 'Y' should be a matrix.") }
  if (base::ncol(X) != base::ncol(Y)){
    stop("* wasserstein : input 'X' and 'Y' should be of same dimension.")
  }
  
  # Number of observation in each matrix
  m <- base::nrow(X)
  n <- base::nrow(Y)
  
  wxname <-  paste0("'",deparse(substitute(wx)),"'")
  wyname <- paste0("'",deparse(substitute(wy)),"'")
  fname  <- "wasserstein"
  
  # Weight normalization
  par_wx <- valid_single_marginal(wx, m, fname)
  par_wy <- valid_single_marginal(wy, n, fname)
  
  # Cost matrix
  dist_mat  <- compute_pdist_simplex_fast(X, Y)
  
  # Solve the optimal transport problem
  wass_lp(dxy = dist_mat, wx = par_wx, wy = par_wy, p = 2)
}