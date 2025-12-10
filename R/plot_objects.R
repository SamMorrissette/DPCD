#' Plot the Object Configuration
#'
#' Generates a plot of the posterior mean of the latent coordinates (`x`) from a DPCD model fit, aligned to a specified target matrix using a Procrustes transformation.
#'
#' @param mcmc_samples An object of class `mcmc` or `mcmc.list` containing posterior samples from a DPCD model fit using [run_dpcd()]. Variable `x` must be included in the output parameters.
#' @param target_matrix A matrix used as the target for aligning the posterior latent coordinates (`x`) via a Procrustes transformation.
#' @param show_clusters Logical argument indicating whether to colour points by their cluster membership. If `TRUE`, then `z` must be monitored in `mcmc_samples`.
#' @param ... Additional arguments to be passed to `plot()` (2 dimensions) or `pairs()` (higher dimensions).
#' @details Since the latent coordinates are non-identifiable due to invariance of Euclidean distances to rotation, reflection, and translation, this function first aligns the posterior samples of `x` to a specified target matrix using a Procrustes transformation. Then, it computes the posterior mean of the aligned latent coordinates and generates a plot. If `show_clusters` is set to `TRUE`, points are coloured according to their cluster memberships, which is estimated through maximizing the posterior expected adjusted Rand index (Fritsch and Ickstadt, 2009).
#'
#' @references Fritsch, Arno & Ickstadt, Katja. (2009). An Improved Criterion for Clustering Based on the Posterior Similarity Matrix. Bayesian Analysis. 4. 10.1214/09-BA414.
#' @returns A scatter plot (for 2-dimensional latent space) or pairs plot (for higher dimensions) of the object configuration.
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(10*2), ncol = 2)
#' dis_matrix <- dist(x)
#' mcmc_samples <- run_dpcd("UU", dis_matrix, p = 2, niter = 10000, nburn = 2000)
#' target_matrix <- cmdscale(dis_matrix, k = 2)
#' plot_x(mcmc_samples, target_matrix, show_clusters = TRUE)
#' }
#' @importFrom graphics pairs
#' @export
plot_objects <- function(mcmc_samples, target_matrix, show_clusters = TRUE, ...) {

  if (!inherits(mcmc_samples, "mcmc.list") & !inherits(mcmc_samples, "mcmc")) {
    stop("`mcmc_samples` must be an object of class 'mcmc.list' or 'mcmc'.")
  }

  if (inherits(mcmc_samples, "mcmc.list")) {
    full_samples <- do.call(rbind, mcmc_samples)
  } else {
    full_samples <- mcmc_samples
  }

  x_samples <- full_samples[, startsWith(colnames(full_samples), "x")]

  if (length(x_samples) == 0) {
    stop("Variable `x` must be monitored in `mcmc_samples`." )
  }

  z_samples <- full_samples[, startsWith(colnames(full_samples), "z")]

  if (show_clusters & length(z_samples) == 0) {
    stop("Variable `z` must be monitored in `mcmc_samples` when show_clusters = TRUE." )
  }

  x_cols <- colnames(x_samples)
  dims <- sub(".*,(\\s*[0-9]+\\s*)\\]", "\\1", x_cols)
  dims <- as.integer(dims)
  latent_dim <- max(dims)

  if (ncol(target_matrix) != latent_dim) {
    stop("`target_matrix` and `x` have different dimensions.")
  }

  iters <- nrow(full_samples)

  x_transformed <- matrix(NA, nrow = iters, ncol = ncol(x_samples))
  for (i in 1:iters) {
    x_mat <- matrix(x_samples[i,], ncol = latent_dim, byrow = FALSE)
    x_transformed[i,] <- c(procrustes(target_matrix, x_mat))
  }
  final_mat <- matrix(colMeans(x_transformed), ncol = latent_dim, byrow = FALSE)
  colnames(final_mat) <- paste0("Dim", 1:latent_dim)

  if (show_clusters == TRUE) {
    cl <- extract_clusters(mcmc_samples)
    if (latent_dim == 2) {
      plot(final_mat, col = cl, pch = 19, ...)
    } else if (latent_dim > 2) {
      pairs(final_mat, col = cl, pch = 19, ...)
    } else {
      stop("Latent dimension must be at least 2.")
    }
  }

  if (show_clusters == FALSE) {
    if (latent_dim == 2) {
      plot(final_mat, pch = 19, ...)
    } else if (latent_dim > 2) {
      pairs(final_mat, pch = 19, ...)
    } else {
      stop("Latent dimension must be at least 2.")
    }
  }
}
