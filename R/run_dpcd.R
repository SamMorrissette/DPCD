#' Run Dirichlet Process Clustering with Dissimilarities
#'
#' @description This function fits an infinite mixture model to dissimilarity data using a Dirichlet Process prior. The model is constructed and MCMC sampling is performed using the [`nimble`](https://r-nimble.org/) package.
#'
#' @param dis_matrix A distance structure such as that returned by [stats::dist] or a full symmetric matrix containing the dissimilarities.
#' @param p The dimension of the space in which the objects are embedded. Must be at least 2.
#' @param hyper_params A named list of hyperparameter values. See details for more information.
#' @param init_params A named list of initial values for model parameters. See details for more information.
#' @param output_params A character vector of model parameters to save in the output. See details for more information.
#' @param scale Logical indicating whether to scale the distance matrix so that its maximum value is 1.
#' @param nchains Number of MCMC chains to run.
#' @param niter Number of MCMC iterations to run.
#' @param nburn Number of MCMC burn-in iterations.
#' @param trunc_value The truncation level for the stick-breaking representation of the Dirichlet process.
#' @param ... Additional arguments passed to [`nimbleMCMC()`] from the `nimble` package.
#'
#' @details Test
#' @returns Test
#' @examples
#' x <- matrix(rnorm(10*2), ncol = 2)
#' dis_matrix <- dist(x)
#' mcmc_samples <- run_dpcd(dis_matrix, p = 2, niter = 1000, nburn = 200)
#'
#' @importFrom stats cmdscale dist
#' @importFrom utils modifyList
#' @export

run_dpcd <- function(dis_matrix,
                     p,
                     trunc_value = 10,
                     hyper_params,
                     init_params,
                     output_params,
                     scale = TRUE,
                     nchains = 1, niter = 10000, nburn = 0, ...) {

  if (!is.matrix(dis_matrix) & !inherits(dis_matrix, "dist")) {
    stop("`dis_matrix` must be a distance matrix.")
  }

  d_obs <- as.matrix(dis_matrix)

  if (!isSymmetric(d_obs)) {
    stop("`dis_matrix` must be a symmetric matrix.")
  }

  if (any(d_obs < 0)) {
    stop("`dis_matrix` must contain only non-negative elements.")
  }

  if (any(diag(d_obs) != 0)) {
    stop("The diagonal elements of `dis_matrix` must be 0.")
  }

  p <- as.integer(p)

  if (p <= 1) {
    stop("p must be a positive integer greater than 1.")
  }

  default_hyper_params <- list("alpha_0" = 1,
                               "a_0" = 1,
                               "b_0" = 1,
                               "lambda" = 1,
                               "mu_0" = rep(0, p),
                               "Psi_0" = diag(p),
                               "nu_0" = p+2)

  if (missing(hyper_params)) {
    hyper_params <- default_hyper_params
  } else {
    if (!is.list(hyper_params)) {
      stop("`hyper_params` must be a list")
    }

    if (!any(names(hyper_params) %in% names(default_hyper_params))) {
      stop("`hyper_params` must be a named list containing values for one or more hyperparameters.")
    }
    hyper_params <- modifyList(default_hyper_params, hyper_params)
  }

  n <- nrow(d_obs)
  constants <- list("p" = p,
                    "alpha_0" = hyper_params$alpha_0,
                    "a_0" = hyper_params$a_0,
                    "b_0" = hyper_params$b_0,
                    "lambda" = hyper_params$lambda,
                    "Psi_0" = hyper_params$Psi_0,
                    "mu_0" = hyper_params$mu_0,
                    "nu_0" = hyper_params$nu_0,
                    n = n,
                    N = trunc_value)

  mds_init <- cmdscale(d_obs, k = p)
  default_inits <- list(
    x = mds_init,
    delta = as.matrix(dist(mds_init)),
    z = sample(1:constants$N, constants$n, replace = TRUE),
    beta = rep(0.5, constants$N - 1),
    Sigma = array(diag(p), dim = c(p,p,constants$N)),
    mu = matrix(0, nrow = constants$N, ncol = p),
    sigma_sq = 1
  )

  if (missing(init_params)) {
    init_params <- default_inits
  } else {
    if (!is.list(init_params)) {
      stop("`init_params` must be a list")
    }

    if (!any(names(init_params) %in% names(default_inits))) {
      stop("`init_params` must be a named list containing initial values for one or more model parameters.")
    }
    init_params <- modifyList(default_inits, init_params)
  }

  default_output_params <- c("x", "z", "pi", "mu", "Sigma", "sigma_sq")
  if (missing(output_params)) {
    monitors <- default_output_params
  } else {
    if (!is.vector(output_params)) {
      stop("`output_vars` must be a character vector.")
    }

    param_names <- c("beta", "pi", "z", "mu", "Sigma", "sigma_sq", "x", "delta")
    if (!all(output_params %in% param_names)) {
      stop("`output_vars` must be a character vector containing valid model parameter names.")
    }
    monitors <- output_params
  }

  if (scale == TRUE) {
    scalar <- 1 / max(d_obs)
    d_obs <- scalar * d_obs
  }

  data <- list(d_obs = d_obs)

  run_uu(data, constants, init_params,
         nchains, niter, nburn,
         monitors,
         ...)
}
