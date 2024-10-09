#' Calculation of the beta-binomial joint density function
#'
#' This function returns a functor that calculates the joint density of two
#' random variables following a beta-binomial distribution.
#'
#' @param x1 Numeric vector. Observations for the first sample.
#' @param x2 Numeric vector. Observations for the second sample.
#' @param n1 Integer. Number of trials in the first sample.
#' @param n2 Integer. Number of trials in the second sample.
#' @param u Numeric. Upper argument for the generalized hypergeometric function.
#' @param l Numeric. Lower argument for the generalized hypergeometric function.
#' @param a0 Numeric. First shape parameter of the beta-binomial distribution.
#' @param a1 Numeric. Second shape parameter of the beta-binomial distribution.
#' @param a2 Numeric. Third shape parameter of the beta-binomial distribution.
#' @param i Integer. Index for the observation to calculate in both samples. Default is 1.
#' @param equal_thetas Logical. If TRUE, theta1 and theta2 will be forced to be equal. Default is FALSE.
#' @param use_log Logical. If TRUE, the function returns the log of the density. Default is FALSE.
#'
#' @return A function that computes the beta-binomial joint density for the given parameters.
#'
#' @examples
#' dens_fun <- densBB_functor(x1 = c(0, 1), x2 = c(1, 1), n1 = 10, n2 = 12,
#'                            u = 1, l = 1, a0 = 0.5, a1 = 1, a2 = 1, i = 1)
#' dens_fun(0.6, 0.7)
#'
#' @importFrom hypergeo genhypergeo
#' @export


densBB_functor <- function(x1, x2, n1, n2, u, l, a0, a1, a2, i = 1, equal_thetas = FALSE, use_log = FALSE) {

  a <- sum(c(a0,a1,a2))

  densBB <- function(theta1, theta2 = NULL) {
    if (equal_thetas) {
      theta2 <- theta1
    }

    if (is.null(theta2)) {
      stop("theta2 debe ser proporcionado cuando equal_thetas es FALSE")
    }

    log_f <- lgamma(n1 + a) + lgamma(n2 + a) -
      lgamma(x1[i] + a1) - lgamma(n1 - x1[i] + a - a1) -
      lgamma(x2[i] + a2) - lgamma(n2 - x2[i] + a - a2) -
      log(genhypergeo(U = u, L = l, check_mod = TRUE, z = 1)) +
      (a1 + x1[i] - 1) * log(theta1) + (a2 + a0 + (n1 - x1[i]) - 1) * log(1 - theta1) +
      (a2 + x2[i] - 1) * log(theta2) + (a1 + a0 + (n2 - x2[i]) - 1) * log(1 - theta2) -
      (a) * log(1 - theta1 * theta2)

    if (!use_log) {
      return(exp(log_f))
    } else {
      return(log_f)
    }
  }
  return(densBB)
}
