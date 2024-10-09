#' Compute adaptive cut-off value k
#'
#' This function computes the value k that minimizes the sum of the errors
#'
#'
#' @param k Numeric. The parameter value for which the approximation is computed.
#' @param x1 Numeric vector. Observations for the first sample.
#' @param x2 Numeric vector. Observations for the second sample.
#' @param n1 Integer. Number of trials for the first sample.
#' @param n2 Integer. Number of trials for the second sample.
#' @param alpha0 Numeric. First shape parameter.
#' @param alpha1 Numeric. Second shape parameter.
#' @param alpha2 Numeric. Third shape parameter.
#' @param j Integer. Index for the iteration loop.
#' @param N1 Integer. Number of iterations for the first loop.
#' @param N2 Integer. Number of iterations for the second loop.
#' @param n_bootstrap Integer. Number of bootstrap samples.
#' @param seed Integer. Seed for random number generation. Default is `NULL`, which generates a random seed.
#'
#' @return Numeric. The Monte Carlo approximation, or stops with an error message if something goes wrong.
#'
#'
#' @examples
#' alf_integ2(k = 0.5, x1 = c(1, 2), x2 = c(2, 3), n1 = 10, n2 = 12, alpha0 = 1, alpha1 = 2, alpha2 = 3, j = 1, N1 = 100, N2 = 100, n_bootstrap = 1000)
#'
#' @importFrom hypergeo genhypergeo
#' @export




select_k <- function(x1, x2, n1, n2, alpha0, alpha1, alpha2, i, j, N2, n_bootstrap, k = seq(0.0000, 1.0000, length.out = 100), seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  alpha_500 <- full_evidence(x1, x2, n1, n2, alpha0, alpha1, alpha2, i, j, N2, n_bootstrap, k, seed)

  f1 <- spline(k, alpha_500, n = 1000, method = "fmm", xmin = min(k), xmax = max(k))
  alpha_500_s <- smooth.spline(f1$x, f1$y, df = 6)$y

  beta_500 <- full_evidence_beta(x1, x2, n1, n2, alpha0, alpha1, alpha2, i, j, N2, n_bootstrap, k, seed)

  f2 <- spline(k, beta_500, n = 1000, method = "fmm", xmin = min(k), xmax = max(k))
  beta_500_s <- smooth.spline(f2$x, f2$y, df = 5)$y

  sumerrors_s <- alpha_500_s + beta_500_s
  k1 <- f1$x
  kop_s <- min(k1[sumerrors_s == min(sumerrors_s)])

  sumerrors <- alpha_500 + beta_500
  kop <- min(k[sumerrors == min(sumerrors)])

  return(kop)
}
