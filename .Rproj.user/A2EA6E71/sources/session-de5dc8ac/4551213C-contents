#' Calculate Evidence for Multiple Grids of Alpha
#'
#' This function calculates the evidence (`alpha_results_all_grids`) for multiple grids of values `k`
#' using parallel computation. It uses the `furrr` and `future` packages to distribute the computation
#' across multiple cores. The function calls `alf_integ2` for each value in the grid, passing the required parameters.
#'
#' @param x1 Numeric vector. Observations for the first sample.
#' @param x2 Numeric vector. Observations for the second sample.
#' @param n1 Integer. Number of trials for the first sample.
#' @param n2 Integer. Number of trials for the second sample.
#' @param alpha0 Numeric. First shape parameter.
#' @param alpha1 Numeric. Second shape parameter.
#' @param alpha2 Numeric. Third shape parameter.
#' @param N1 Integer. Number of iterations for the first loop.
#' @param N2 Integer. Number of iterations for the second loop.
#' @param n_bootstrap Integer. Number of bootstrap samples.
#' @param seed Integer. Seed for random number generation. Default is `NULL`, which generates a random seed.
#'
#' @return A list where each element corresponds to the results of evidence calculations for a specific grid of `k`.
#'
#' @examples
#' full_evidence_alpha(x1 = c(1, 2), x2 = c(2, 3), n1 = 10, n2 = 12, alpha0 = 1, alpha1 = 2, alpha2 = 3, N1 = 100, N2 = 100, n_bootstrap = 1000)
#'
#' @importFrom hypergeo genhypergeo
#'
#' @export
full_evidence_alpha <- function(x1, x2, n1, n2, alpha0, alpha1, alpha2, N1, N2, n_bootstrap, seed = NULL) {

  library(furrr)
  library(future)

  m_grids_k <- list(
    round(seq(0.0000, 0.0909, length.out = 10), 4),
    round(seq(0.1010, 0.1919, length.out = 10), 4),
    round(seq(0.2020, 0.2929, length.out = 10), 4),
    round(seq(0.3030, 0.3939, length.out = 10), 4),
    round(seq(0.4040, 0.4949, length.out = 10), 4),
    round(seq(0.5051, 0.5960, length.out = 10), 4),
    round(seq(0.6061, 0.6970, length.out = 10), 4),
    round(seq(0.7071, 0.7980, length.out = 10), 4),
    round(seq(0.8081, 0.8990, length.out = 10), 4),
    round(seq(0.9091, 1.0000, length.out = 10), 4)
  )

  j <- seq(1, N2)

  # Planificación paralela
  cores <- 10
  plan(multisession, workers = cores)

  # Función para calcular por grilla
  calculate_for_grid <- function(m_grid_k) {
    future_map_dbl(
      m_grid_k,
      ~mean(sapply(j, function(i) mc_integral_alpha(k = .x, x1 = x1, x2 = x2, n1 = n1, n2 = n2,
                                             alpha0 = alpha0, alpha1 = alpha1, alpha2 = alpha2, j = i,
                                             N1 = N1, N2 = N2, n_bootstrap = n_bootstrap, seed = seed)),
            na.rm = TRUE),
      .options = furrr_options(seed = TRUE)
    )
  }

  # Calcular los resultados para todas las grillas
  alpha_results_all_grids <- lapply(m_grids_k, calculate_for_grid)

  return(alpha_results_all_grids)
}
