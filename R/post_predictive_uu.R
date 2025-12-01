post_predictive_uu <- function(dis_matrix,
                               p,
                               mcmc_samples,
                               nsim = 1000,
                               scale = TRUE,
                               plot = TRUE) {
  # Insert validity check to make sure that all of the input variables are present...

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

  if (scale == TRUE) {
    scalar <- 1 / max(d_obs)
    d_obs <- scalar * d_obs
  }

  n <- nrow(d_obs)
  m <- n * (n-1) / 2
  x_samples <- mcmc_samples[, startsWith(colnames(mcmc_samples), "x")]
  sigma_samples <- sqrt(mcmc_samples[, "sigma_sq"])

  d_pred <- matrix(data = NA, nrow = nsim, ncol = m)

  if (nrow(x_samples) >= nsim) {
    sample_idx <- sample(1:nrow(x_samples), nsim, replace = FALSE)
  } else {
    stop("There are not enough posterior draws in `mcmc_samples`. Please choose a lower `nsim` value.")
  }


  x_sub_samples <- x_samples[sample_idx, , drop = FALSE]
  sigma_sub_samples <- sigma_samples[sample_idx]
  for (i in 1:nsim) {
    delta_vec <- c(dist(matrix(x_sub_samples[i, ], nrow = n, ncol = p)))
    d_pred[i, ] <- truncnorm::rtruncnorm(m,
                                         mean = delta_vec,
                                         sd = sigma_sub_samples[i],
                                         a = 0)
  }

  if (plot == TRUE) {
    d_obs_vec <- c(d_obs[lower.tri(d_obs)])
    p <- bayesplot::ppc_dens_overlay(d_obs_vec, d_pred) +
      ggtitle("Posterior predictive check")
    print(p)
  }

  return(d_pred)
}
