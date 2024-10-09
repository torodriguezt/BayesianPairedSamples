#' Calculate Error Based on Evidence
#'
#' This function calculates an error metric (`error2`) based on the evidence calculated by the `evidence` function under type II error.
#' It includes checks for `NA` values in the calculated evidence (`ev10`) and verification value (`f_ver`),
#' and throws an error if any of them are `NA`.
#'
#' @param x1 Numeric vector. Observations for the first sample.
#' @param x2 Numeric vector. Observations for the second sample.
#' @param n1 Integer. Number of trials for the first sample.
#' @param n2 Integer. Number of trials for the second sample.
#' @param alpha0 Numeric. First shape parameter.
#' @param alpha1 Numeric. Second shape parameter.
#' @param alpha2 Numeric. Third shape parameter.
#' @param i Integer. Index for the observation in both samples.
#' @param j Integer. Index for Theta generation.
#' @param N2 Integer. Number of samples for `Theta1`.
#' @param n_bootstrap Integer. Number of bootstrap samples for variational Bayes.
#' @param k Numeric. Threshold for calculating the indicator `Ind`.
#' @param seed Integer. Seed for random number generation. Default is `NULL`, which generates a random seed.
#'
#' @return Numeric. The calculated error `error2`.
#'
#' @examples
#' error_one(x1 = c(1, 2), x2 = c(2, 3), n1 = 10, n2 = 12, alpha0 = 1, alpha1 = 2, alpha2 = 3, i = 1, j = 1, N2 = 100, n_bootstrap = 1000, k = 0.5)
#'
#' @importFrom hypergeo genhypergeo
#' @importFrom MCMCpack MCMCmetrop1R
#' @export



error_two <- function(x1, x2, n1, n2, alpha0, alpha1, alpha2, i, j, N2, n_bootstrap, k, seed = NULL) {


  ##No estoy muy seguro de como cambiar esta parte

  priorB2<- function(theta) {

    if( any(theta > 1) | any(theta < 0) ) {
      return(-Inf)
    }

    a<-alpha0+alpha1+alpha2

    log_f<-lgamma(a)-
      lgamma(alpha0)-lgamma(alpha1)-
      lgamma(alpha2)+
      (alpha1-1)*log(theta[1])+(alpha2+alpha1-1)*log(1-theta[1])+
      (alpha2-1)*log(theta[2])+(alpha1+alpha1-1)*log(1-theta[2])-
      (a)*log(1-theta[1]*theta[2])
    return(log_f)
  }

  theta.samp <- MCMCmetrop1R(priorB2, theta.init=c(0.3,0.3),
                             thin=10, mcmc=200000, burnin=5000,
                             # tune=c(1, 1),
                             verbose=500, logfun=TRUE,
                             force.samp=T,
                             optim.lower=c(0.01, 0.01),
                             optim.upper = c(0.99, 0.99),
                             optim.method = "L-BFGS-B")

  Theta11<-theta.samp[, 1][seq(20,20000,40)]
  Theta22<-theta.samp[, 2][seq(20,20000,40)]
  N2=length(Theta11)
  N1=length(Theta22)

  ev10 <- evidence_beta(x1, x2, n1, n2, alpha0, alpha1, alpha2, i, j, N2, n_bootstrap, theta1 = Theta11, theta2 = Theta22 ,seed)

  Ind <- I(ev10 > k)


  error2 <- f_ver * Ind

  return(error2)
}



