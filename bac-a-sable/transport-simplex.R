#' Gaussian-based transport of a single observation from group 0 to group 1
#' on \eqn{S_d}, using additive log ratio transform
#'
#' @param x Value to transport.
#' @param n_interp Number of points for the interpolation.
#' @param A Symmetric positive matrix from OT.
#' @param m0 Average of transformed values in group 0.
#' @param m1 Average of transformed values in group 1.
#'
#' @importFrom compositions alr alrInv
#' @export
#' @md
transport_x_alr <- function(x,
                            n_interp = 31,
                            A,
                            m0,
                            m1) {
  
  z <- as.numeric(alr(as.numeric(x)))
  vec_t <- seq(0, 1, length = n_interp)
  transp_x <- matrix(NA, n_interp, length(x))
  for(i in 1:n_interp) {
    t <- vec_t[i]
    transp_z <- (1 - t) * z + t * (m1 + A %*% (z - m0))
    transp_z = as.numeric(transp_z)
    transp_x[i,] = alrInv(transp_z)
  }
  transp_x
}

#' Gaussian-based transport of a single observation from group 0 to group 1
#' on \eqn{S_d}, using centered log ratio transform
#'
#' @param x Value to transport.
#' @param n_interp Number of points for the interpolation.
#' @param A Symmetric positive matrix from OT.
#' @param m0 Average of transformed values in group 0.
#' @param m1 Average of transformed values in group 1.
#'
#' @importFrom compositions clr clrInv
#' @export
#' @md
transport_x_clr <- function(x,
                            n_interp = 31,
                            A,
                            m0,
                            m1) {
  
  z <- as.numeric(clr(as.numeric(x)))[1:(length(x)-1)]
  vec_t <- seq(0, 1, length = n_interp)
  transp_x <- matrix(NA, n_interp, length(x))
  for(i in 1:n_interp) {
    t <- vec_t[i]
    transp_z <- (1 - t) * z + t * (m1 + A %*% (z - m0))
    transp_z = as.numeric(transp_z)
    transp_z = c(transp_z, -sum(transp_z))
    transp_x[i,] = clrInv(transp_z)
  }
  transp_x
}

#' Gaussian-based transport of a single observation from group 0 to group 1
#' on \eqn{S_d}, using isometric log ratio transform
#'
#' @param x Value to transport.
#' @param n_interp Number of points for the interpolation.
#' @param A Symmetric positive matrix from OT.
#' @param m0 Average of transformed values in group 0.
#' @param m1 Average of transformed values in group 1.
#'
#' @importFrom compositions ilr ilrInv
#' @export
#' @md
transport_x_ilr <- function(x,
                            n_interp = 31,
                            A,
                            m0,
                            m1) {
  
  z <- as.numeric(ilr(as.numeric(x)))
  vec_t <- seq(0, 1, length = n_interp)
  transp_x <- matrix(NA, n_interp, length(x))
  for(i in 1:n_interp) {
    t <- vec_t[i]
    transp_z <- (1 - t) * z + t * (m1 + A %*% (z - m0))
    transp_z = as.numeric(transp_z)
    transp_x[i,] = ilrInv(transp_z)
  }
  transp_x
}

#' Learn OT mapping from group 0 to group 1 based on the representation of the
#' composition data from the simplex (\eqn{S_d}) to a \eqn{d-1} Euclidean space.
#'
#' @param X0 Data frame with observations in group 0.
#' @param X1 Data frame with observations in group 1.
#' @param isomorphism Isomorphism used to map the composition data from the
#'  simplex (\eqn{S_d}) to a \eqn{d-1} Euclidean space. Three possibilities: `"alr"`
#'  (additive log ratio transform), `"clr"` (centered log ratio transform),
#'  and `"ilr"` (isometric log ratio transform).
#'
#'  @returns A list with the following elements:
#'  * `m0`, `m1`: empirical mean of \eqn{z_0} and \eqn{z_1} (vectors of orthonormal
#'    coordinates).
#'  * `S0`, `S1`: empirical covariance matrices of \eqn{z_0} and \eqn{z_1}.
#'  * `A`: symmetric positive matrix that satisfies
#'    \eqn{\bf{A}\bf{\Sigma}_{0}\bf{A}=\bf{\Sigma}_{1}}
#'  * `isomorphism`: name of the isomorphism used (alr, clr, irl).
#'
#' @importFrom compositions alr clr ilr alrInv clrInv ilrInv
#' @importFrom stats var
#' @importFrom expm sqrtm
#' @export
#' @md
#' @examples
#' # First three columns: probabilities of being of class A, B, or C.
#' # Last column: group (0 or 1)
#' data(toydataset)
#' X0 <- toydataset[toydataset$group == 0, c("A", "B", "C")]
#' X1 <- toydataset[toydataset$group == 1, c("A", "B", "C")]
#' ot_mapping <- get_ot_mapping(X0 = X0, X1 = X1, isomorphism = "clr")
#' ot_mapping$A
get_ot_mapping <- function(X0,
                           X1,
                           isomorphism = c("clr", "alr", "ilr")) {
  
  isomorphism <- match.arg(isomorphism)
  
  if (isomorphism == "alr") {
    Z0 <- alr(X0)
    Z1 <- alr(X1)
    transport_f_x <- transport_x_alr
  } else if (isomorphism == "clr") {
    Z0 <- matrix(clr(X0), ncol = ncol(X0))
    Z1 <- matrix(clr(X1), ncol = ncol(X1))
    Z0 <- Z0[, 1:(ncol(Z0)-1)]
    Z1 <- Z1[, 1:(ncol(Z1)-1)]
    transport_f_x <- transport_x_clr
  } else {
    Z0 <- ilr(X0)
    Z1 <- ilr(X1)
    transport_f_x <- transport_x_ilr
  }
  
  # empirical mean in each group
  m0 <- apply(Z0, 2, mean)
  m1 <- apply(Z1, 2, mean)
  # empirical variance in each group
  S0 <- var(Z0)
  S1 <- var(Z1)
  
  A <- solve(sqrtm(S0)) %*%
    sqrtm(sqrtm(S0) %*% S1 %*% (sqrtm(S0))) %*%
    solve(sqrtm(S0))
  
  list(
    m0 = m0,
    m1 = m1,
    S0 = S0,
    S1 = S1,
    A = A,
    isomorphism = isomorphism
  )
  
}

#' Transport compositional data from group 0 to group 1
#'
#' @param X0 Data frame with observations in group 0.
#' @param X1 Data frame with observations in group 1.
#' @param isomorphism Isomorphism used to map the composition data from the
#'  simplex (\eqn{S_d}) to a \eqn{d-1} Euclidean space. Three possibilities: `"alr"`
#'  (additive log ratio transform), `"clr"` (centered log ratio transform),
#'  and `"ilr"` (isometric log ratio transform).
#' @param n_interp Number of steps in the interpolation (default to 1: no
#'  interpolation).
#'
#' @returns A tibble with the transported values. If `n_interp` is larger than
#'  1, the result also contains a list with the interpolated values, in the
#'  `"interpolated"` attribute ; else, the attribute is `NULL`. The attribute
#'  `"ot_mapping"` stores the mapping.
#'
#' @references  McCann, Robert J. 1997. "A Convexity Principle for Interacting
#' Gases." Advances in Mathematics 128 (1): 153–79.
#'
#' @seealso [transport_simplex_new()]
#'
#' @importFrom rlang set_names
#' @importFrom tibble as_tibble
#' @importFrom purrr list_rbind map
#' @importFrom dplyr slice_tail
#' @export
#' @md
#' @examples
#' # First three columns: probabilities of being of class A, B, or C.
#' # Last column: group (0 or 1)
#' data(toydataset)
#' X0 <- toydataset[toydataset$group == 0, c("A", "B", "C")]
#' X1 <- toydataset[toydataset$group == 1, c("A", "B", "C")]
#'
#' # Transport only, from group 0 to group 1, using centered log ratio transform:
#' transp <- transport_simplex(X0 = X0, X1 = X1, isomorphism = "clr")
#' head(transp)
#'
#' # If we want to transport new points:
#' new_obs <- data.frame(A = c(.2, .1), B = c(.6, .5), C = c(.2, .4))
#' transport_simplex_new(transport = transp, newdata = new_obs)
#'
#' # If we want to get interpolated values using McCann (1997) displacement
#' # interpolation: (here, with 5 intermediate points)
#' transp_with_interp <- transport_simplex(
#'   X0 = X0, X1 = X1, isomorphism = "clr", n_interp = 5
#' )
#' interpolated(transp_with_interp)[[1]] # first obs
#' interpolated(transp_with_interp)[[2]] # second obs
#'
#' # And displacement interpolation for the new obs:
#' transp_new_obs_with_interp <- transport_simplex_new(
#'   transport = transp, newdata = new_obs, n_interp = 5
#' )
#' interpolated(transp_new_obs_with_interp)[[1]] # first new obs
#' interpolated(transp_new_obs_with_interp)[[1]] # second new obs
transport_simplex <- function(X0,
                              X1,
                              isomorphism = c("clr", "alr", "ilr"),
                              n_interp = 1) {
  isomorphism <- match.arg(isomorphism)
  ot_mapping <- get_ot_mapping(X0 = X0, X1 = X1, isomorphism = isomorphism)
  # empirical mean in each group
  m0 <- ot_mapping$m0
  m1 <- ot_mapping$m1
  # empirical variance in each group
  S0 <- ot_mapping$S0
  S1 <- ot_mapping$S1
  
  A <- ot_mapping$A
  
  if (isomorphism == "alr") {
    transport_f_x <- transport_x_alr
  } else if (isomorphism == "clr") {
    transport_f_x <- transport_x_clr
  } else {
    transport_f_x <- transport_x_ilr
  }
  
  transported <- map(
    1:nrow(X0),
    ~{
      transp_val <- transport_f_x(
        x = X0[.x, ],  n_interp = n_interp, A = A, m0 = m0, m1 = m1
      )
      colnames(transp_val) <- colnames(X0)
      as_tibble(transp_val)
    }
  )
  
  if (n_interp == 1) {
    transported_val <- transported |> list_rbind()
    interpolated <- NULL
  } else {
    transported_val <- map(transported, ~slice_tail(.x, n = 1)) |> list_rbind()
    interpolated <- transported
  }
  structure(
    transported_val,
    interpolated = interpolated,
    ot_mapping = ot_mapping
  )
}

#' Transport new compositional data from group 0 to group 1 using a previously
#' learned mapping.
#'
#' @param transport Previously learned mapping.
#' @param newdata New data frame with composition data.
#' @param n_interp Number of steps in the interpolation (default to 1: no
#'  interpolation).
#'
#' @returns A tibble with the transported values. If `n_interp` is larger than
#'  1, the result also contains a list with the interpolated values, in the
#'  `"interpolated"` attribute ; else, the attribute is `NULL`. The attribute
#'  `"ot_mapping"` stores the mapping.
#'
#' @references  McCann, Robert J. 1997. "A Convexity Principle for Interacting
#' Gases." Advances in Mathematics 128 (1): 153–79.
#'
#' @seealso [transport_simplex()]
#'
#' @importFrom rlang set_names
#' @importFrom tibble as_tibble
#' @importFrom purrr list_rbind map
#' @importFrom dplyr slice_tail
#' @export
#' @md
#' @examples
#' # First three columns: probabilities of being of class A, B, or C.
#' # Last column: group (0 or 1)
#' data(toydataset)
#' X0 <- toydataset[toydataset$group == 0, c("A", "B", "C")]
#' X1 <- toydataset[toydataset$group == 1, c("A", "B", "C")]
#'
#' # Transport only, from group 0 to group 1, using centered log ratio transform:
#' transp <- transport_simplex(X0 = X0, X1 = X1, isomorphism = "clr")
#' head(transp)
#'
#' # If we want to transport new points:
#' new_obs <- data.frame(A = c(.2, .1), B = c(.6, .5), C = c(.2, .4))
#' transp_new_obs <- transport_simplex_new(transport = transp, newdata = new_obs)
#' transp_new_obs
#'
#' # If we want to get interpolated values using McCann (1997) displacement
#' # interpolation: (here, with 5 intermediate points)
#' transp_new_obs_with_interp <- transport_simplex_new(
#'   transport = transp, newdata = new_obs, n_interp = 5
#' )
#' interpolated(transp_new_obs_with_interp)[[1]] # first new obs
#' interpolated(transp_new_obs_with_interp)[[2]] # second new obs
transport_simplex_new <- function(transport,
                                  newdata,
                                  n_interp = 1) {
  
  ot_mapping <- attr(transport, "ot_mapping")
  
  # empirical mean in each group
  m0 <- ot_mapping$m0
  m1 <- ot_mapping$m1
  # empirical variance in each group
  S0 <- ot_mapping$S0
  S1 <- ot_mapping$S1
  
  A <- ot_mapping$A
  
  isomorphism <- ot_mapping$isomorphism
  
  if (isomorphism == "alr") {
    transport_f_x <- transport_x_alr
  } else if (isomorphism == "clr") {
    transport_f_x <- transport_x_clr
  } else {
    transport_f_x <- transport_x_ilr
  }
  
  transported <- map(
    1:nrow(newdata),
    ~{
      transp_val <- transport_f_x(
        x = newdata[.x, ],  n_interp = n_interp, A = A, m0 = m0, m1 = m1
      )
      colnames(transp_val) <- colnames(newdata)
      as_tibble(transp_val)
    }
  )
  
  if (n_interp == 1) {
    transported_val <- transported |> list_rbind()
    interpolated <- NULL
  } else {
    transported_val <- map(transported, ~slice_tail(.x, n = 1)) |> list_rbind()
    interpolated <- transported
  }
  structure(
    transported_val,
    interpolated = interpolated,
    ot_mapping = ot_mapping
  )
}

#' Extract the interpolated values of transported vectors of compositional data
#'
#' @param x Transported compositional data.
#' @returns A list with the interpolated values. Each element contains a tibble
#'  giving the interpolated values for an observation.
#' @export
#' @seealso [transport_simplex()]
interpolated <- function(x) {
  attr(x, "interpolated")
}


#' Cost function for optimal transport on the unit simplex
#' @noRd
d_s <- function(x, y) {
  d <- length(x)
  log(mean(y / x)) - mean(log(y / x))
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
compute_pdist_simplex <- function(X, Y) {
  M <- matrix(NA, nrow(X), nrow(Y))
  for (i in 1:nrow(X)) {
    for (j in 1:nrow(Y)) {
      M[i, j] <- d_s(X[i, ], Y[j, ])
    }
  }
  M
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
  cxy    <- (dxy)
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
  dist_mat  <- compute_pdist_simplex(X, Y)
  
  # Solve the optimal transport problem
  wass_lp(dist_mat, par_wx, par_wy, p = 2)
}


#' Builds counterfactuals of observations from group 0 to group 1 using an
#' estimated mapping.
#'
#' @param mapping Wasserstein mapping between the two sets of probability
#'  vectors
#' @param X0 Matrix (or data frame) with vectors of probabilities in group 0.
#' @param X1 Matrix (or data frame) with vectors of probabilities in group 1.
#'
#' @returns A matrix with the counterfactual values using the mapping.
#'
#' @export
#'
#' @examples
#' data(toydataset)
#' # Observations in group 0
#' X0 <- toydataset[toydataset$group == 0, c("A", "B", "C")]
#' # Observations in group 1
#' X1 <- toydataset[toydataset$group == 1, c("A", "B", "C")]
#'
#' # Optimal Transport using Linear Programming:
#' mapping <- wasserstein_simplex(as.matrix(X0), as.matrix(X1))
#'
#' # The counterfactuals of observations of group 0 in group 1
#' counterfactuals_0_1 <- counterfactual_w(mapping, X0, X1)
#' head(counterfactuals_0_1)
counterfactual_w <- function(mapping,
                             X0,
                             X1) {
  if (!is.matrix(X0)) X0 <- as.matrix(X0)
  if (!is.matrix(X1)) X1 <- as.matrix(X1)
  M0 <- mapping$plan * nrow(X0)
  
  M0 %*% X1
}
