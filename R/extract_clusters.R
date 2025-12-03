#' Title
#'
#' @param mcmc_samples
#'
#' @returns
#' @export
#'
#' @examples
extract_clusters <- function(mcmc_samples) {

  if (!inherits(mcmc_samples, "mcmc.list") & !inherits(mcmc_samples, "mcmc")) {
    stop("`mcmc_samples` must be an object of class 'mcmc.list' or 'mcmc'.")
  }

  if (inherits(mcmc_samples, "mcmc.list")) {
    full_samples <- do.call(rbind, mcmc_samples)
  } else if (inherits(mcmc_samples, "mcmc")) {
    full_samples <- mcmc_samples
  }

  z_samples <- full_samples[, startsWith(colnames(full_samples), "z")]

  if (length(z_samples) == 0) {
    stop("Variable `z` must be monitored in `mcmc_samples`." )
  }

  psm <- mcclust::comp.psm(z_samples)
  cluster_labels <- mcclust::maxpear(psm)

  cluster_labels$cl
}
