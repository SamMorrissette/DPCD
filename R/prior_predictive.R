#' Prior Predictive Check
#'
#' @inheritParams run_dpcd
#' @param nsim Number of datasets to simulate from the prior predictive distribution.
#' @param plot Logical indicating whether to plot the simulated datasets. See details for more information.
#'
#' @returns
#' @export
#'
#' @examples
prior_predictive <- function(model_name = c("UU", "EU", "UD", "ED", "US", "ES"),
                             dis_matrix,
                             p,
                             trunc_value = 10,
                             hyper_params = NULL,
                             init_params = NULL,
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
    scalar <- 1 / max(d_obs)
    d_obs <- scalar * d_obs
  }

  model_config <- get_config(model_name, d_obs, p, trunc_value, hyper_params, init_params)

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
