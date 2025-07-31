#' Dataset Generation with a Sensitive Attribute
#'
#' @description
#' This function generates a dataset where covariates \eqn{X = (X1, X2)} are drawn from bivariate Gaussian distributions that depend on a binary sensitive attribute \eqn{S}. The binary outcome \eqn{Y} is generated based on \eqn{X} and \eqn{S} through a logistic model.
#'
#' @details
#' For \eqn{S = 0}, \eqn{X} is sampled from a multivariate normal distribution with mean vector \eqn{[-1, -1]} and covariance matrix \deqn{\begin{bmatrix} 1.2^2 & 1.2^2 \cdot 0.5 \\ 1.2^2 \cdot 0.5 & 1.2^2 \end{bmatrix}.}
#'
#' For \eqn{S = 1}, \eqn{X} is sampled from a multivariate normal distribution with mean vector \eqn{[1.5, 1.5]} and covariance matrix \deqn{\begin{bmatrix} 0.9^2 & 0.9^2 \cdot -0.4 \\ 0.9^2 \cdot -0.4 & 0.9^2 \end{bmatrix}.}
#'
#' The logistic model uses a linear predictor \eqn{\eta} to determine the probability \eqn{p} of \eqn{Y = 1}:
#' * For \eqn{S = 0}, \eqn{\eta = 0.6 \cdot X1 + 0.2 \cdot X2}.
#' * For \eqn{S = 1}, \eqn{\eta = 0.4 \cdot X1 + 0.6 \cdot X2}.
#'
#' \eqn{Y} is sampled from a Bernoulli distribution with probability \eqn{p}, where \deqn{p = \frac{\exp(\eta)}{1 + \exp(\eta)}.}
#'
#' The final dataset includes \eqn{S}, \eqn{X1}, \eqn{X2}, and \eqn{Y}.
#'
#' @param n_0 Number of observations in group S=0 (Default to 100).
#' @param n_1 Number of observations in group S=1 (Default to 100).
#' @param mu_0 Mean of the bivariate Gaussian in group S=0. (Default to -1).
#' @param mu_1 Mean of the bivariate Gaussian in group S=1. (Default to 1.5).
#' @param sigma_0 A symmetric positive-definite matrix representing the
#'        variance-covariance matrix of the distribuion in group S=0.
#' @param sigma_1 A symmetric positive-definite matrix representing the
#'        variance-covariance matrix of the distribuion in group S=1.
#'
#' @importFrom mnormt rmnorm
#' @importFrom stats rbinom
#' @importFrom tibble as_tibble
#'
#' @examples
#' sim_dat <- simul_dataset()
#' head(sim_dat)
#' @export
#' @md
simul_dataset <- function(n_0 = 100,
                          n_1 = 100,
                          mu_0 = rep(-1, 2),
                          mu_1 = rep(1.5, 2),
                          sigma_0 = matrix(c(1, .5, .5, 1) * 1.2^2, 2, 2),
                          sigma_1 = matrix(c(1, -.4, -.4, 1) * .9^2, 2, 2)) {
  # Covariates from bivariate Gaussians
  X0 <- mnormt::rmnorm(n_0, mu_0, sigma_0)
  X1 <- mnormt::rmnorm(n_1, mu_1, sigma_1)
  D_SXY_0 <- data.frame(S = 0, X1 = X0[, 1], X2 = X0[, 2])
  D_SXY_1 <- data.frame(S = 1, X1 = X1[, 1], X2 = X1[, 2])
  # Outcome: depends on X1 and X2. First, predictor
  eta_0 <- (D_SXY_0$X1 * 1.2 + D_SXY_0$X2 / 2 * .8) / 2
  eta_1 <- (D_SXY_1$X1 * .8 + D_SXY_1$X2 / 2 * 1.2) / 2
  # Apply logistic function to get values in [0,1]
  p_0 <- exp(eta_0) / (1 + exp(eta_0))
  p_1 <- exp(eta_1) / (1 + exp(eta_1))
  # Random draw in Binomial distribution
  D_SXY_0$Y <- rbinom(n_0, size = 1, prob = p_0)
  D_SXY_1$Y <- rbinom(n_1, size = 1, prob = p_1)
  D_SXY <- rbind(D_SXY_0, D_SXY_1)
  as_tibble(D_SXY)
}
