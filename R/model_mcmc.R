#' MCMC Sampling for a Specified Nimble Model
#'
#' An internal function to run the unequal unrestricted (UU) model. This function is called by [run_dpcd()].
#'
#' @param data A symmetric dissimilarity matrix.
#' @param constants A named list of model constant values (including hyperparameters).
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
model_mcmc <- function(model_name,
                       data,
                       constants,
                       init_params,
                       WAIC,
                       nchains,
                       niter,
                       nburn,
                       output_params,
                       ...) {

  model_code <- switch(
    model_name,
    "UU" = unequal_unrestricted_model,
    "EU" = equal_unrestricted_model,
    "UD" = unequal_diagonal_model,
    "ED" = equal_diagonal_model,
    "US" = unequal_spherical_model,
    "ES" = equal_spherical_model,
  )
  N <- p <- lambda <- Sigma <- n <- x <- NULL
  nimble::nimbleMCMC(code = model_code,
                     constants = constants,
                     data = data,
                     inits = init_params,
                     monitors = output_params,
                     WAIC = WAIC,
                     nchains = nchains,
                     niter = niter,
                     nburnin = nburn,
                     samplesAsCodaMCMC = TRUE,
                     ...)
}




