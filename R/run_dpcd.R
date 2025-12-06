#' Run Dirichlet Process Clustering with Dissimilarities
#'
#' @description This function fits an infinite mixture model to dissimilarity data using a Dirichlet Process prior. The model is constructed and MCMC sampling is performed using the [`nimble`](https://r-nimble.org/) package. Currently, there are six different models available.
#'
#' @param model_name The DPCD model to fit. Must be one of "UU" (unequal unrestricted), "EU" (equal unrestricted), "UD" (unequal diagonal), "ED" (equal diagonal), "US" (unequal spherical), or "ES" (equal spherical). See details for a brief description of each model.
#' @param dis_matrix A distance structure such as that returned by [stats::dist] or a full symmetric matrix containing the dissimilarities.
#' @param p The dimension of the space in which the objects are embedded. Must be at least 2.
#' @param hyper_params A named list of hyperparameter values. See details for more information.
#' @param init_params A named list of initial values for model parameters. See details for more information.
#' @param output_params A character vector of model parameters to save in the output. See details for more information.
#' @param scale Logical argument indicating whether to scale the dissimilarities so that the maximum value is 1.
#' @param nchains Number of MCMC chains to run.
#' @param niter Number of MCMC iterations to run.
#' @param nburn Number of MCMC burn-in iterations.
#' @param trunc_value The truncation level for the stick-breaking representation of the Dirichlet process.
#' @param WAIC Logical argument indicating whether to compute the Watanabe-Akaike Information Criterion (WAIC) for model comparison.
#' @param ... Additional arguments passed to [`nimbleMCMC()`] from the `nimble` package.
#'
#' @details
#' Dirichlet Process Clustering with Dissimilarities (DPCD) models dissimilarity
#' data using an infinite mixture model with a Dirichlet Process prior. The six
#' available covariance structures for mixture components are:
#'
#' - **"UU"**: Unequal Unrestricted — each component has its own unrestricted covariance matrix.
#' - **"EU"**: Equal Unrestricted — components share a common unrestricted covariance matrix.
#' - **"UD"**: Unequal Diagonal — each component has its own diagonal covariance matrix.
#' - **"ED"**: Equal Diagonal — components share a common diagonal covariance matrix.
#' - **"US"**: Unequal Spherical — each component has its own spherical covariance matrix.
#' - **"ES"**: Equal Spherical — components share a common spherical covariance matrix.
#'
#' The `hyper_params` list allows users to specify custom hyperparameter values.
#' Some hyperparameters are common across all models, while others depend on the
#' selected covariance structure.
#'
#' **Common hyperparameters:**
#' - `alpha_0`: Concentration parameter for the Dirichlet Process prior.
#' - `a_0`, `b_0`: Shape and scale parameters for the Inverse-Gamma prior on the
#'   measurement error parameter.
#' - `lambda`: Scaling parameter for the prior on component means.
#' - `mu_0`: Mean vector for the prior on component means.
#'
#' **Model-specific hyperparameters:**
#' - `nu_0` and `Psi_0`(degrees of freedom and scale matrix for the Inverse-Wishart prior) - UU and EU only.
#' - `alpha_tau` and `beta_tau` (shape and scale
#'   parameters for the Inverse-Gamma prior) - UD, ED, US, and ES only.
#'
#' The `init_params` list allows users to supply initial values for model
#' parameters to assist MCMC convergence. The following parameters may be
#' initialized:
#' - `x`: `n × p` matrix of latent positions.
#' - `sigma_sq`: Scalar measurement error variance.
#' - `mu`: `trunc_value × p` matrix of component means.
#' - `Sigma`: `p × p` covariance matrix.
#' - `tau_sq`: Scalar variance parameter (for "US" and "ES" only).
#' - `tau_vec`: Length-`p` variance vector (for "UD" and "ED" only).
#' - `beta`: Length `trunc_value-1` vector of stick-breaking weights.
#' - `z`: Length-`n` vector of cluster assignments.
#'
#' Default values are used for both `hyper_params` and `init_params` if none are
#' supplied.
#'
#' The `output_params` vector specifies which model parameters should be saved in
#' the MCMC output. Valid names include `"beta"`, `"pi"`, `"z"`, `"mu"`,
#' `"Sigma"`, `"sigma_sq"`, `"x"`, and `"delta"`.
#'
#' @returns Posterior samples are returned a coda `mcmc` object, unless `nchains > 1`, in which case the posterior samples are returend as a coda `mcmc.list` object. If `WAIC = TRUE`, a named list is returned containing the posterior samples and the WAIC value.
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(10*2), ncol = 2)
#'
#' # Fit the unequal unrestricted model with default settings
#' dis_matrix <- dist(x)
#' mcmc_samples <- run_dpcd("UU", dis_matrix, p = 2, niter = 10000, nburn = 2000)
#' summary(mcmc_samples)
#'
#' # Fit the equal spherical model with custom hyperparameters and initial values
#' custom_hyper_params <- list(alpha_tau = 0.01, beta_tau = 0.01)
#' custom_init_params <- list(sigma_sq = 0.5)
#' mcmc_samples_es <- run_dpcd("ES", dis_matrix, p = 2,
#'                             hyper_params = custom_hyper_params,
#'                             init_params = custom_init_params,
#'                             niter = 10000, nburn = 2000, WAIC = TRUE)
#'}
#'
#' @importFrom stats cmdscale dist
#' @importFrom utils modifyList
#' @export

run_dpcd <- function(model_name = c("UU", "EU", "UD", "ED", "US", "ES"),
                     dis_matrix,
                     p,
                     trunc_value = 15,
                     hyper_params = NULL,
                     init_params = NULL,
                     output_params = c("x", "z", "pi", "mu", "Sigma", "sigma_sq"),
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

  param_names <- c("beta", "pi", "z", "mu", "Sigma", "sigma_sq", "x", "delta")
  if (!all(output_params %in% param_names)) {
    stop("`output_vars` must be a character vector containing valid model parameter names.")
  }

  if (scale == TRUE) {
    d_obs <- d_obs / max(d_obs)
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
