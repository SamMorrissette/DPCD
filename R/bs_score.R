#' Calculate the Bayesian Silhouette Score
#'
#' @description This function calculates the Bayesian Silhouette (BS) Score for a DPCD model fit using posterior MCMC samples. The BS score can be used to evaluate the clustering quality of a fit and to compare different models.
#' @param mcmc_samples An object of class `mcmc` or `mcmc.list` containing posterior samples from a DPCD model fit using [run_dpcd()]. Variables `x` and `z` must be included in the output parameters.
#'
#' @details The Bayesian Silhouette Score is computed by calculating the silhouette score for each MCMC iteration based on the latent positions (`x`) and cluster assignments (`z`). The silhouette score measures how similar an object is to its own cluster compared to other clusters. The BS score is then obtained by averaging the silhouette scores across all MCMC iterations. Higher values of the BS score indicate a higher-quality DPCD model in terms of its clustering structure.
#' @returns A numeric value representing the average silhouette score across all MCMC iterations.
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(10*2), ncol = 2)
#' dis_matrix <- dist(x)
#' mcmc_samples <- run_dpcd("UU", dis_matrix, p = 2, niter = 10000, nburn = 2000)
#' bs_score(mcmc_samples)
#' }
#'
#' @importFrom cluster silhouette
#' @export

bs_score <- function(mcmc_samples) {
  if (!inherits(mcmc_samples, "mcmc.list") & !inherits(mcmc_samples, "mcmc")) {
    stop("`mcmc_samples` must be an object of class 'mcmc.list' or 'mcmc'.")
  }

  if (inherits(mcmc_samples, "mcmc.list")) {
    full_samples <- do.call(rbind, mcmc_samples)
  } else {
    full_samples <- mcmc_samples
  }

  x_samples <- full_samples[, startsWith(colnames(full_samples), "x")]
  z_samples <- full_samples[, startsWith(colnames(full_samples), "z")]

  if (length(x_samples) == 0) {
    stop("Variable `x` must be monitored in `mcmc_samples`." )
  }

  if (length(z_samples) == 0) {
    stop("Variable `z` must be monitored in `mcmc_samples`.")
  }

  iter <- nrow(z_samples)
  x_cols <- colnames(x_samples)
  dims <- sub(".*,(\\s*[0-9]+\\s*)\\]", "\\1", x_cols)
  dims <- as.integer(dims)
  latent_dim <- max(dims)

  bs <- rep(NA, iter)
  for (t in 1:iter) {
    x_matrix <- matrix(x_samples[t, ], ncol = latent_dim, byrow = FALSE)
    z <- z_samples[t, ]
    dist_mat <- dist(x_matrix)
    num_clusters <- length(unique(z))
    if (num_clusters > 1) {
      sil <- cluster::silhouette(z, dist_mat)
      bs[t] <- mean(sil[,'sil_width'])
    } else {
      bs[t] <- NA
    }
  }
  if (any(is.na(bs))) {
    warning(paste(sum(is.na(bs)), "silhouette values were NA. This is usually caused by iterations where all objects were placed into one cluster."))
  }
  mean(bs, na.rm = TRUE)
}
