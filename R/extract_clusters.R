#' Extract clusters from MCMC samples
#'
#' This function extracts estimated cluster memberships from MCMC samples obtained from a DPCD model fit.
#'
#' @param mcmc_samples An object of class `mcmc` or `mcmc.list` containing posterior samples from a DPCD model fit using [run_dpcd()]. The variable `z` must be included in the output parameters.
#'
#' @details
#' This function uses the cluster membership variable, `z`, from the provided MCMC samples to compute the posterior similarity matrix (PSM) based on the sampled cluster assignments. Using the PSM, it then determines the estimated cluster memberships by maximizing the posterior expected adjusted Rand index, following the method of Fritsch and Ickstadt (2009).
#'
#' @references Fritsch, Arno & Ickstadt, Katja. (2009). An Improved Criterion for Clustering Based on the Posterior Similarity Matrix. Bayesian Analysis. 4. <doi:10.1214/09-BA414>.
#'
#' @seealso [mcclust::maxpear()]
#' @returns A vector of labels that indicate the estimated cluster membership for each observation.
#'
#' @examples
#' extract_clusters(mcmc_example)
#'
#'
#' @importFrom mcclust comp.psm maxpear
#' @export

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
