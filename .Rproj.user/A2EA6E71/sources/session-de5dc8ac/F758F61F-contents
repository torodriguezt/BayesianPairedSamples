#' Calculate Evidence for the Beta-Binomial Model
#'
#' This function calculates the evidence for a beta-binomial model using variational inference under error type I
#' via Stan and genetic algorithms to optimize parameters. It returns the calculated evidence (`ev10`).
#'
#' @param x1 Numeric vector. Observations for the first sample.
#' @param x2 Numeric vector. Observations for the second sample.
#' @param n1 Integer. Number of trials for the first sample.
#' @param n2 Integer. Number of trials for the second sample.
#' @param alpha0 Numeric. First shape parameter.
#' @param alpha1 Numeric. Second shape parameter.
#' @param alpha2 Numeric. Third shape parameter.
#' @param i Integer. Index for the observation in both samples.
#' @param j Integer. Index for Theta sampling.
#' @param N2 Integer. Number of samples for `Theta1`.
#' @param n_bootstrap Integer. Number of bootstrap samples for variational Bayes.
#' @param seed Integer. Seed for random number generation. Default is `NULL`, which generates a random seed.
#'
#' @return Numeric. The calculated evidence `ev10`.
#'
#' @examples
#' evidence(x1 = c(1, 0), x2 = c(1, 1), n1 = 2, n2 = 2, alpha0 = 0.6, alpha1 = 0.65, alpha2 = 0.7, i = 1, j = 1, N2 = 100, n_bootstrap = 1000)
#'
#' @importFrom hypergeo genhypergeo
#' @importFrom rstan extract
#'
#' @export
#'
#'
#'
evidence <- function(x1, x2, n1, n2, alpha0, alpha1, alpha2, i, j, N2, n_bootstrap, seed = NULL) {

  if (is.null(seed)) {
    seed <- sample.int(.Machine$integer.max, 1)
  }

  f2 <- function(X1, X2, theta1) {
    exp(dbinom(X1, n1, theta1, log = TRUE) + dbinom(X2, n2, theta1, log = TRUE))
  }

  alpha <- sum(c(alpha0, alpha1, alpha2))
  Theta1 <- rbeta(N2, 1, 1)

  f_ver <- f2(x1[i], x2[i], Theta1[j])

  u <- c(alpha, x1[i] + alpha1, x2[i] + alpha2)
  l <- c(n1 + alpha, n2 + alpha)

  stan_model <- rstan::stan_model(file = 'inst/stan/BBpost3.stan', verbose = TRUE)

  GA_ev <- ga(type = "real-valued",
              fitness = function(x) {
                densBB_functor(x1, x2, n1, n2, u, l, alpha0, alpha1, alpha2, i, equal_thetas = TRUE)(x[1])
              },
              lower = c(0), upper = c(0.99),
              popSize = 50, maxiter = 1000, run = 100,
              monitor = FALSE)

  supremo_ga_ev <- GA_ev@solution[1, ] |> as.numeric()

  n <- as.integer(c(n1, n2))
  X <- c(x1[i], x2[i])
  P <- length(X)

  stan_data <- list(P = P, X = X, n = n, alpha1 = alpha0, alpha2 = alpha1, alpha3 = alpha2)

  sink(nullfile())
  stan_vb <- rstan::vb(object = stan_model, data = stan_data, seed = seed,
                       output_samples = n_bootstrap)
  sink()

  Thetas <- cbind(
    extract(stan_vb)$Theta[, 1] |> as.vector(),
    extract(stan_vb)$Theta[, 2] |> as.vector()
  )

  densBB <- densBB_functor(x1, x2, n1, n2, u, l, alpha0, alpha1, alpha2, i, equal_thetas = TRUE)

  ev10 <- mean(apply(Thetas, 1, function(t) {
    I(densBB(t[1], t[2]) < densBB(supremo_ga_ev))
  }), na.rm = TRUE)

  return(ev10)
}
