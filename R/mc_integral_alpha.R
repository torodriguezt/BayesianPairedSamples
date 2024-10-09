#' Compute Monte Carlo Approximation for Evidence
#'
#' This function computes a Monte Carlo approximation for evidence by repeatedly calling `error_one`
#' over a sequence of iterations. If any errors occur during the computation, it stops with a descriptive error message.
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
#' @examples
#' alf_integ2(k = 0.5, x1 = c(1, 2), x2 = c(2, 3), n1 = 10, n2 = 12, alpha0 = 1, alpha1 = 2, alpha2 = 3, j = 1, N1 = 100, N2 = 100, n_bootstrap = 1000)
#'
#' @importFrom hypergeo genhypergeo
#' @export
mc_integral_alpha <- function(k, x1, x2, n1, n2, alpha0, alpha1, alpha2, j, N1, N2, n_bootstrap, seed = NULL) {

  i_seq <- seq(1, N1)

  MC1 <- tryCatch({
    mean(sapply(i_seq, function(i) error_one(x1 = x1, x2 = x2, n1 = n1, n2 = n2, alpha0 = alpha0, alpha1 = alpha1,
                                             alpha2 = alpha2, i = i, j = j, N2 = N2, n_bootstrap = n_bootstrap, k = k, seed = seed)), na.rm = TRUE)
  }, error = function(e) {
    stop("Error en el cálculo de Monte Carlo: ", e$message)
  })

  return(MC1)
}
