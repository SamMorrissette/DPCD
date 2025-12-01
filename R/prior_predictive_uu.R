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
prior_predictive_uu <- function(dis_matrix,
                                p,
                                trunc_value = 10,
                                hyper_params,
                                init_params,
                                output_params,
                                scale = TRUE,
                                nsim = 1000,
                                plot = TRUE) {

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

  if (scale == TRUE) {
    scalar <- 1 / max(d_obs)
    d_obs <- scalar * d_obs
  }

  empty_data <- lapply(list(d_obs = d_obs), function(x) array(NA, dim = dim(x)))

  pp_samples <- nimbleMCMC(
    code = unequal_unrestricted_model,
    constants = constants,
    data = empty_data,
    inits = init_params,
    monitors = c("x", "sigma_sq"),
    nchains = 1,
    niter = nsim,
    nburnin = 0,
    samplesAsCodaMCMC = TRUE
  )

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
      ggplot2::ggtitle("Prior predictive check")
    print(p)
  }

  return(d_pred)
}
