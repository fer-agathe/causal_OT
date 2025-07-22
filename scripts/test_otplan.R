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
  wass_lp_sinkhorn(dxy = dist_mat, wx = par_wx, wy = par_wy, p = 2)
}

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

wass_lp_sinkhorn <- function(dxy, wx, wy, p = 2) {
  stopifnot(all(abs(sum(wx) - 1) < 1e-8), all(abs(sum(wy) - 1) < 1e-8))
  
  # Compute transport plan via Sinkhorn algorithm
  gamma <- T4transport::sinkhornD(dxy, p = 2, wx = wx, wy = wy, lambda = 0.1)
  
  list(distance = gamma$distance, plan = gamma$plan)
}

wass_lp_fast <- function(dxy, wx, wy, p = 2) {
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
