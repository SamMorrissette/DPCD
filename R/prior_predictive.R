#' Prior Predictive Check
#'
#' This function simulates dissimilarities from the prior predictive distribution of a specified DPCD model and optionally plots the density of the simulated dissimilarities against the observed dissimilarities.
#'
#' @param model_name The DPCD model from which to draw prior predictive samples. Must be one of "UU", "EU", "UD", "ED", "US", or "ES".
#' @param nsim Number of datasets to simulate from the prior predictive distribution.
#' @param plot Logical argument indicating whether to plot the simulated dissimilarities against the observed dissimilarities. See details for more information.
#' @inheritParams run_dpcd
#'
#' @details A prior predictive check is used to assess if datasets drawn from the prior predictive distribution are consistent with the observed data. Most of the mass of the prior predictive distribution should be placed on plausible values of the dissimilarities, while little or no mass should be placed on implausible values.
#'
#' If `plot = TRUE`, a plot is created to compare the density of the observed dissimilarities to the densities of the dissimilarities simulated from the prior predictive distribution using `bayesplot::ppc_dens_overlay()`.
#'
#' See [run_dpcd()] for details on the DPCD models and hyperparameters.
#' @seealso [run_dpcd()]
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A. (2019).
#' Visualization in Bayesian workflow. Journal of the Royal Statistical Society A,
#' 182(2), 389–402. https://doi.org/10.1111/rssa.12378
#'
#' @returns A matrix of simulated dissimilarities from the prior predictive distribution with `nsim` rows and `n * (n-1) / 2` columns, where `n` is the number of objects (i.e. the number of rows/columns of `dis_matrix`).
#' @export
#'
#' @examples
#' \donttest{
#' ppc <- prior_predictive(dis_mat_example, "UU", p = 2, nsim = 100, plot = TRUE)
#' }
#'
#' @import ggplot2
#' @importFrom bayesplot ppc_dens_overlay
#' @importFrom truncnorm rtruncnorm
#'
prior_predictive <- function(dis_matrix,
                             model_name = c("UU", "EU", "UD", "ED", "US", "ES"),
                             p = 2,
                             trunc_value = 15,
                             hyper_params = NULL,
                             scale = TRUE,
                             nsim = 1000,
                             plot = TRUE) {

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

  if (scale == TRUE) {
    d_obs <- d_obs / max(d_obs)
  }

  model_config <- get_config(model_name, d_obs, p, trunc_value, hyper_params, user_inits = NULL)

  empty_data <- list(d_obs = array(NA, dim = dim(d_obs)))

  pp_samples <- model_mcmc(
    model_name,
    empty_data,
    model_config$constants,
    model_config$init_params,
    WAIC = FALSE,
    nchains = 1,
    niter = nsim,
    nburn = 0,
    output_params = c("x", "sigma_sq"))

  n <- nrow(d_obs)
  m <- n * (n-1) / 2
  x_samples <- pp_samples[, startsWith(colnames(pp_samples), "x")]
  sigma_samples <- sqrt(pp_samples[, "sigma_sq"])

  d_pred <- matrix(data = NA, nrow = nsim, ncol = m)
  for (i in 1:nsim) {
    delta_vec <- c(dist(matrix(x_samples[i, ], nrow = n, ncol = p)))
    d_pred[i, ] <- truncnorm::rtruncnorm(m,
                                         mean = delta_vec,
                                         sd = sigma_samples[i],
                                         a = 0)
  }

  if (plot == TRUE) {
    d_obs_vec <- c(d_obs[lower.tri(d_obs)])
    p <- bayesplot::ppc_dens_overlay(d_obs_vec, d_pred) +
      ggplot2::ggtitle(paste("Prior predictive check for model", model_name))
    print(p)
  }

  d_pred
}
