# remotes::install_github(repo = "fer-agathe/transport-simplex")
library(transportsimplex)

# Data----

set.seed(1234)
n <- 1000
X1 <- rnorm(n)
X2star <- X1 + rnorm(n)
X2 <- cut(X2star, c(-Inf, -1, +1, +Inf) * .7, c("A", "B", "C"))
X3 <- X1 + (X2 == "B") + rnorm(n)
eta <- X1 + (X2 == "A") - 0.4
probaT <- exp(eta) / (1 + exp(eta))
T <- rbinom(n, size = 1, prob = probaT)
Y <- 1 + 2 * X1 + rnorm(n) + T * (X3 + (X2 == "C"))
base <- data.frame(
  X1 = X1,
  X2 = as.factor(X2),
  X3 = X3,
  T  = T,
  Y = Y
)

# Propensities for the categorical variable----
# using a Multinomial Logistic Regression

library(nnet)
library(splines)
propensities_reg <- multinom(X2 ~ X1 + X3 + bs(X1) + bs(X3), data = base)
prob2 <- predict(propensities_reg,newdata=base,type="probs")
base$probA <- prob2[,1]
base$probB <- prob2[,2]
base$probC <- prob2[,3]

# Split data by treatment
base0 = base[base$T == 0,]
base1 = base[base$T == 1,]

# Compare means
apply(base1[,c("X1","X3","T","Y")],2,mean)
apply(base0[,c("X1","X3","T","Y")],2,mean)

# Start counterfactual construction for treated
base1_counterfactual <- base1
base1_counterfactual$T <- 0

# Counterf. X1----
# Construct OT map for X1 (numeric variable)
F1_1 <- Vectorize(function(x) mean(base1$X1 <= x))
Q1_0 <- Vectorize(function(u) as.numeric(quantile(base0$X1, u)))
Transp1 <- function(x)  Q1_0( F1_1(x))
base1_counterfactual$X1 <- Transp1(base1$X1)

# Counterf. X2----

# Use optimal transport on the simplex to map class probabilities
# We loop over each treated to build a **local** counterfactual.
cli::cli_progress_bar("X2: treated --> untreated", total = nrow(base1))
for (i in 1:nrow(base1)) {
  h <- .5
  x1 <- base1$X1[i]
  Tx1 <- base1_counterfactual$X1[i]
  
  # Kernel weights in the neighborhood
  poids1 <- dnorm(base1$X1, x1, h)
  poids0 <- dnorm(base0$X1, Tx1, h)
  poids1 <- poids1 / sum(poids1)
  poids0 <- poids0 / sum(poids0)
  # Identify 50 closest neighbors
  # r1 = which(rank(poids1) > nrow(base1) - 50)
  # r0 = which(rank(poids0) > nrow(base0) - 50)
  r1 <- order(poids1, decreasing = TRUE)[1:50]
  r0 <- order(poids0, decreasing = TRUE)[1:50]
  ir1 <- which(r1 == i)
  
  # Optimal Transport using Linear Programming:
  # Wasserstein distance between the set of probability vectors in
  # first group (treated) and second group (untreated)
  Wi = wasserstein_simplex(
    X = as.matrix(base1[r1,c("probA","probB","probC")]),
    Y = as.matrix(base0[r0,c("probA","probB","probC")]),
    wx = poids1[r1],
    wy = poids0[r0]
  )
  
  # Most likely match under the transport plan
  Pi <- base0[r0, c("probA","probB","probC")][which.max(Wi$plan[ir1, ]), ]
  
  base1_counterfactual$probA[i] <- Pi[["probA"]]
  base1_counterfactual$probB[i] <- Pi[["probB"]]
  base1_counterfactual$probC[i] <- Pi[["probC"]]
  cli::cli_progress_update()
}

# OT to assign new X2 values to counterfactuals
library(transport)
proportions0 <- table(base0$X2) / nrow(base0) 
proportions1 <- table(base1$X2) / nrow(base1)

samples <- as.data.frame(base1_counterfactual[, c("probA","probB","probC")])
# Target probabilities for each class of X2
p_target <- proportions0
p_target

# n <- nrow(samples)
# colnames(samples) <- c("A", "B", "C")
# rownames(samples) <- 1:n
# vertices <- data.frame(
#   A = c(1, 0, 0),
#   B = c(0, 1, 0),
#   C = c(0, 0, 1)
# )
# cost_matrix <- as.matrix(dist(rbind(samples, vertices))^2)
# cost_matrix <- cost_matrix[1:n, (n+1):(n+3)]
# colnames(cost_matrix) <- c("A", "B", "C")
# mass_source <- rep(1/n, n)
# mass_target <- as.numeric(p_target)
# names(mass_target) <- c("A", "B", "C")
# plan <- transport::transport(a = mass_source, 
#                              b = mass_target, 
#                              costm = cost_matrix, 
#                              method = "shortsimplex")
# assignments <- rep(NA, n)
# for (i in 1:nrow(plan)) {
#   assignments[plan$from[i]] <- plan$to[i]
# }
# samples$category <- factor(assignments, labels = c("A", "B", "C"))

#' OT for categorical variable, from source distribution to target 
#' probabilities.
#' 
#' @param x Data frame with individuals (rows) from the source distribution,
#'  where the columns give the probabilities associated with each class.
#' @param labels Vector of labels of the categorical variable.
#' @param p Vector of target probabilities. If omitted, uniform weights are 
#'  used.
#' 
get_assignment <- function(x,
                           labels = NULL,
                           p = NULL) {
  
  n_labels <- ncol(x)
  n <- nrow(x)
  if (is.null(p)) p <- rep(1, n_labels) / n_labels # Uniform weights
  
  # Unit vectors
  vertices <- diag(n_labels)
  colnames(vertices) <- colnames(x)
  # source weights
  mass_source <- rep(1 / n, n)
  # target weights
  mass_target <- as.numeric(p)
  
  # Cost matrix (squared Euclidean distance)
  cost_matrix <- as.matrix(dist(rbind(x, vertices))^2)
  cost_matrix <- cost_matrix[1:n, (n + 1):(n + n_labels)]
  
  # Assign eah observation to one vertex
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
  assignments <- max.col(mass_matrix, ties.method = "first")
  
  factor(assignments, labels = labels)
}

base1_counterfactual$X2 <- get_assignment(
  x = samples, labels = c("A", "B", "C"), p = p_target
)

# Counterf. X3----

library(DescTools)

for (i in 1:nrow(base1)) {
  h <- .5
  x1 <- base1$X1[i]
  x2 <- base1$X2[i]
  Tx1 <- base1_counterfactual$X1[i]  
  Tx2 <- base1_counterfactual$X2[i]
  
  # Weights in the neighborhood
  poids1 <- dnorm(base1$X1[base1$X2 == x2], x1, h)
  poids0 <- dnorm(base0$X1[base0$X2 == Tx2], Tx1, h)
  poids1 <- poids1 / sum(poids1)
  poids0 <- poids0 / sum(poids0)
  
  # Weighted quantile transport
  F3_1 <- Vectorize(
    function(x) weighted.mean(base1$X1[base1$X2 == x2] <= x, weights = poids1)
  )
  Q3_0 <- Vectorize(
    function(u) as.numeric(
      DescTools::Quantile(base0$X1[base0$X2 == Tx2], weights = poids0, u)
    )
  )
  
  Transp3 <- function(x) Q3_0(F3_1(x))
  base1_counterfactual$X3[i] <- Transp3(base1$X3[i])
}

# Compare empirical distributions
table(base0$X2) / nrow(base0) * 100
table(base1_counterfactual$X2) / nrow(base1_counterfactual)*100

# Before/after transport for X1
plot(density(base0$X1), col = "blue")
lines(density(base1$X1), col = "red")
lines(density(base1_counterfactual$X1), col = "purple", lwd = 2, lty = 2)

# Before/after transport for X3 | X2, X1*
plot(density(base0$X3), col = "blue")
lines(density(base1$X3), col = "red")
lines(density(base1_counterfactual$X3), col = "purple", lwd = 2, lty = 2)


# Outcome models
mu_hat_0 = lm(Y~X1+X2+X3+bs(X1)+bs(X3),data=base0)
mu_hat_1 = lm(Y~X1+X2+X3+bs(X1)+bs(X3),data=base1)

# Treatment effects using counterfactuals
TT = predict(mu_hat_1, newdata = base1) - 
  predict(mu_hat_0, newdata = base1_counterfactual)

# CATE vs X1
vx1 <- seq(-1,1,length = 101)
h1 <- function(x1) {
  h <- .5
  poids1 <- dnorm(base1$X1,x1,h)
  poids1 <- poids1 / sum(poids1)
  weighted.mean(TT, poids1)
}
vy1 <- Vectorize(h1)(vx1)
plot(vx1, vy1, type = "l", lwd = 3)
plot(vx1[-1], diff(vy1) / diff(vx1), type = "l", lwd = 3, ylim = c(.5, 1.5))
abline(h = 1)

# CATE vs X3
vx3 <- seq(-1, 1, length = 101)
h3 <- function(x3) {
  h <- .5
  poids1 <- dnorm(base1$X3,x3,h)
  poids1 <- poids1/sum(poids1)
  weighted.mean(TT,poids1)
}
vy3 <- Vectorize(h3)(vx3)
plot(vx3, vy3, type = "l", lwd = 3)
plot(vx3[-1], diff(vy3) / diff(vx3), type = "l", lwd = 3, ylim = c(.5, 1.5))
abline(h = 1)

# ATE by group (X2)
barplot(tapply(TT, base1$X2, mean))


# Brouillon----

# 
# d_s <- function(x, y) {
#   d <- length(x)
#   log(mean(y / x)) - mean(log(y / x))
# }
# compute_pdist_simplex <- function(X, Y) {
#   M <- matrix(NA, nrow(X), nrow(Y))
#   for (i in 1:nrow(X)) {
#     for (j in 1:nrow(Y)) {
#       M[i, j] <- d_s(X[i, ], Y[j, ])
#     }
#   }
#   M
# }
# valid_single_marginal <- function(mvec, M, fname) {
#   dname <- paste0("'", deparse(substitute(mvec)), "'")
#   if ((length(mvec) == 0) && is.null(mvec)) {
#     return(rep(1 / M, M))
#   } else {
#     mvec <- as.vector(mvec)
#     if ((length(mvec) != M) || (any(mvec < 0))) {
#       stop(
#         paste0(
#           "* ", fname, " : ", dname,
#           " should be a nonnegative vector of length ",M,"."
#         )
#       )
#     }
#     return(mvec / base::sum(mvec))
#   }
# }
# wass_lp <- function(dxy,
#                     wx,
#                     wy,
#                     p) {
#   cxy    <- (dxy)
#   m      <- length(wx)
#   ww_m   <- matrix(wx, ncol = 1)
#   n      <- length(wy)
#   ww_n   <- matrix(wy, nrow = 1)
#   ones_m <- matrix(rep(1, n), ncol = 1)
#   ones_n <- matrix(rep(1, m), nrow = 1)
#   plan   <- CVXR::Variable(m, n)
#   
#   wd.obj    <- CVXR::Minimize(CVXR::matrix_trace(t(cxy) %*% plan))
#   wd.const1 <- list(plan >= 0)
#   wd.const2 <- list(plan %*% ones_m == ww_m, ones_n %*% plan == ww_n)
#   wd.prob   <- CVXR::Problem(wd.obj, c(wd.const1, wd.const2))
#   wd.solve  <- CVXR::solve(wd.prob, solver = "OSQP")
#   
#   if (all(wd.solve$status=="optimal")) {
#     # successful
#     gamma <- wd.solve$getValue(plan)
#     value <- (base::sum(gamma * cxy))
#   } else {
#     # failed : use lpsolve
#     cxy <- (dxy)
#     m   <- nrow(cxy)
#     n   <- ncol(cxy)
#     
#     c  <- as.vector(cxy)
#     A1 <- base::kronecker(matrix(1, nrow = 1, ncol = n), diag(m))
#     A2 <- base::kronecker(diag(n), matrix(1, nrow = 1, ncol = m))
#     A  <- rbind(A1, A2)
#     
#     f.obj <- c
#     f.con <- A
#     f.dir <- rep("==", nrow(A))
#     f.rhs <- c(rep(1 / m, m), rep(1 / n, n))
#     f.sol <- (lpSolve::lp("min", f.obj, f.con, f.dir, f.rhs))
#     
#     gamma <- matrix(f.sol$solution, nrow = m)
#     value <- (sum(gamma*cxy)^(1 / p))
#   }
#   list(distance = value, plan = gamma)
# }
# wasserstein_simplex <- function(X,
#                                 Y,
#                                 wx = NULL,
#                                 wy = NULL) {
#   ## CHECK INPUTS
#   if (is.vector(X)) {
#     X <- matrix(X, ncol = 1)
#   }
#   if (is.vector(Y)) {
#     Y <- matrix(Y, ncol = 1)
#   }
#   if (!is.matrix(X)) { stop("* wasserstein : input 'X' should be a matrix.") }
#   if (!is.matrix(Y)) { stop("* wasserstein : input 'Y' should be a matrix.") }
#   if (base::ncol(X) != base::ncol(Y)){
#     stop("* wasserstein : input 'X' and 'Y' should be of same dimension.")
#   }
#   
#   # Number of observation in each matrix
#   m <- base::nrow(X)
#   n <- base::nrow(Y)
#   
#   wxname <-  paste0("'",deparse(substitute(wx)),"'")
#   wyname <-  paste0("'",deparse(substitute(wy)),"'")
#   fname  <- "wasserstein"
#   
#   # Weight normalization
#   par_wx <- valid_single_marginal(wx, m, fname)
#   par_wy <- valid_single_marginal(wy, n, fname)
#   
#   # Cost matrix
#   dist_mat  <- compute_pdist_simplex(X, Y)
#   
#   # Solve the optimal transport problem
#   wass_lp(dist_mat, par_wx, par_wy, p = 2)
# }
# 
