#' Posterior Predictive Check
#'
#' This function simulates dissimilarities from the posterior predictive distribution of a specified DPCD model and optionally plots the density of the simulated dissimilarities against the observed dissimilarities.
#'
#' @param mcmc_samples An object of class `mcmc` or `mcmc.list` containing posterior samples from a DPCD model fit using [run_dpcd()]. Both the latent positions `x` and the error variance `sigma_sq` must be included in `mcmc_samples`.
#' @param nsim Number of datasets to simulate from the posterior predictive distribution.
#' @inheritParams prior_predictive
#' @details A posterior predictive check is used to assess if datasets drawn from the posterior predictive distribution are consistent with the observed data. Posterior predictive checks differ from prior predictive checks in that they incorporate information from the observed data. If the model fits the data well, the observed dissimilarities should look similar to dissimilarities simulated from the posterior predictive distribution.
#'
#' If `plot = TRUE`, a plot is created to compare the density of the observed dissimilarities to the densities of the dissimilarities simulated from the posterior predictive distribution using `bayesplot::ppc_dens_overlay()`.
#'
#' See [run_dpcd()] for details on the DPCD models and hyperparameters.
#' @seealso [run_dpcd()]
#' @references
#' Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., & Gelman, A. (2019).
#' Visualization in Bayesian workflow. Journal of the Royal Statistical Society A,
#' 182(2), 389–402. https://doi.org/10.1111/rssa.12378
#' @returns A matrix of simulated dissimilarities from the posterior predictive distribution with `nsim` rows and `n * (n-1) / 2` columns, where `n` is the number of objects (i.e. the number of rows/columns of `dis_matrix`).
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(10*2), ncol = 2)
#' dis_matrix <- dist(x)
#'
#' # Fit the UU model to simulated data
#' mcmc_samples <- run_dpcd("UU", dis_matrix, p = 2, niter = 10000, nburn = 2000)
#'
#' # Perform a posterior predictive check.
#' ppc <- post_predictive(dis_matrix, mcmc_samples, nsim = 1000, plot = TRUE)
#' }
#' @import ggplot2
#' @importFrom bayesplot ppc_dens_overlay
#' @export
post_predictive <- function(mcmc_samples,
                            dis_matrix,
                            nsim = 1000,
                            scale = TRUE,
                            plot = TRUE) {

  if (!inherits(mcmc_samples, "mcmc.list") & !inherits(mcmc_samples, "mcmc")) {
    stop("`mcmc_samples` must be an object of class 'mcmc.list' or 'mcmc'.")
  }

  if (inherits(mcmc_samples, "mcmc.list")) {
    full_samples <- do.call(rbind, mcmc_samples)
  } else {
    full_samples <- mcmc_samples
  }

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

  x_samples <- full_samples[, startsWith(colnames(full_samples), "x")]
  sigma_samples <- sqrt(full_samples[, "sigma_sq"])

  if (length(x_samples) == 0 || length(sigma_samples) == 0) {
    stop("Variables `x` and `sigma_sq` must both be monitored in `mcmc_samples`.")
  }

  if (scale == TRUE) {
    d_obs <- d_obs / max(d_obs)
  }

  n <- nrow(d_obs)
  m <- n * (n-1) / 2

  x_cols <- colnames(x_samples)
  dims <- sub(".*,(\\s*[0-9]+\\s*)\\]", "\\1", x_cols)
  dims <- as.integer(dims)
  latent_dim <- max(dims)

  d_pred <- matrix(data = NA, nrow = nsim, ncol = m)

  if (nrow(x_samples) >= nsim) {
    sample_idx <- sample(1:nrow(x_samples), nsim, replace = FALSE)
  } else {
    stop("There are not enough posterior draws in `mcmc_samples`. Please choose a lower `nsim` value.")
  }


  x_sub_samples <- x_samples[sample_idx, , drop = FALSE]
  sigma_sub_samples <- sigma_samples[sample_idx]
  for (i in 1:nsim) {
    delta_vec <- c(dist(matrix(x_sub_samples[i, ], nrow = n, ncol = latent_dim)))
    d_pred[i, ] <- truncnorm::rtruncnorm(m,
                                         mean = delta_vec,
                                         sd = sigma_sub_samples[i],
                                         a = 0)
  }

  if (plot == TRUE) {
    d_obs_vec <- c(d_obs[lower.tri(d_obs)])
    p <- bayesplot::ppc_dens_overlay(d_obs_vec, d_pred) +
      ggplot2::ggtitle("Posterior predictive check")
    print(p)
  }

  d_pred
}
