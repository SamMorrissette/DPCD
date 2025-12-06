#' MCMC Sampling for a Specified Nimble Model
#'
#' An internal function to run the unequal unrestricted (UU) model. This is not intended to be called by users. This function is called by [run_dpcd()].
#'
#' @import nimble
#'
#' @keywords internal
#' @noRd
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




