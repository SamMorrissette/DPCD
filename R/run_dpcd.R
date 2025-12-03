#' Run Dirichlet Process Clustering with Dissimilarities
#'
#' @description This function fits an infinite mixture model to dissimilarity data using a Dirichlet Process prior. The model is constructed and MCMC sampling is performed using the [`nimble`](https://r-nimble.org/) package. Currently, there are six different models available.
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

run_dpcd <- function(model_name = c("UU", "EU", "UD", "ED", "US", "ES"),
                     dis_matrix,
                     p,
                     trunc_value = 10,
                     hyper_params = NULL,
                     init_params = NULL,
                     output_params = c("x", "z", "pi", "mu", "Sigma", "sigma_sq", "beta"),
                     scale = TRUE,
                     WAIC = TRUE,
                     nchains = 1, niter = 10000, nburn = 0, ...) {

  model_name <- match.arg(model_name)

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

  param_names <- c("beta", "pi", "z", "mu", "Sigma", "sigma_sq", "x", "delta", "tau", "tau_vec")
  if (!all(output_params %in% param_names)) {
    stop("`output_vars` must be a character vector containing valid model parameter names.")
  }

  if (scale == TRUE) {
    scalar <- 1 / max(d_obs)
    d_obs <- scalar * d_obs
  }

  model_config <- get_config(model_name, d_obs, p, trunc_value, hyper_params, init_params)

  data <- list(d_obs = d_obs)

  model_mcmc(
    model_name,
    data,
    model_config$constants,
    model_config$init_params,
    WAIC,
    nchains,
    niter,
    nburn,
    output_params,
    ...)
}
