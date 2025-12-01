#' Run the Unequal Unrestricted Model
#'
#' An internal function to run the unequal unrestricted (UU) model. This function is called by [run_dpcd()].
#'
#' @param data A symmetric dissimilarity matrix.
#' @param hyper_params A named list of hyperparameter values.
#' @param init_params A named list of initial values for model parameters.
#' @param monitors A character vector of model parameters to save in the output.
#' @param nchains The number of MCMC chains to run.
#' @param niter The number of MCMC iterations to run.
#' @param nburn The number of MCMC burn-in iterations.
#' @param ... Additional arguments passed to `nimbleMCMC()`.
#'
#' @returns MCMC samples as a `coda` object.
#'
#' @import nimble
#'
#' @keywords internal
run_uu <- function(data, hyper_params, init_params,
                                 nchains, niter, nburn,
                                 monitors,
                                 ...) {

  N <- p <- lambda <- Sigma <- n <- x <- NULL
  nimble::nimbleMCMC(code = unequal_unrestricted_model,
                     constants = hyper_params,
                     data = data,
                     inits = init_params,
                     monitors = monitors,
                     nchains = nchains,
                     niter = niter,
                     nburnin = nburn,
                     samplesAsCodaMCMC = TRUE,
                     ...)
}




