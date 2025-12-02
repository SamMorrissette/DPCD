plot_x <- function(mcmc_samples, align_matrix, show_clusters = TRUE) {

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
    stop("Variable `z` must be monitored in `mcmc_samples` when show_clusters is TRUE." )
  }

  x_cols <- colnames(x_samples)
  dims <- sub(".*,(\\s*[0-9]+\\s*)\\]", "\\1", x_cols)
  dims <- as.integer(dims)
  latent_dim <- max(dims)

  if (ncol(align_matrix) != latent_dim) {
    stop("`align_matrix` and `x` have different dimensions.")
  }

  iters <- nrow(full_samples)

  x_transformed <- matrix(NA, nrow = iters, ncol = ncol(x_samples))
  for (i in 1:iters) {
    x_mat <- matrix(x_samples[i,], ncol = latent_dim, byrow = FALSE)
    x_transformed[i,] <- c(procrustes(align_matrix, x_mat))
  }
  final_mat <- matrix(colMeans(x_transformed), ncol = latent_dim, byrow = FALSE)

  if (show_clusters == TRUE) {
    psm <- mcclust::comp.psm(z_samples)
    cluster_labels <- mcclust::maxpear(psm)
    if (latent_dim == 2) {
      plot(final_mat, col = cluster_labels$cl, pch = 19)
    } else if (latent_dim > 2) {
      pairs(final_mat, col = cluster_labels$cl, pch = 19)
    } else {
      stop("Latent dimension must be at least 2.")
    }
  }

  if (show_clusters == FALSE) {
    if (latent_dim == 2) {
      plot(final_mat, pch = 19)
    } else if (latent_dim > 2) {
      pairs(final_mat, pch = 19)
    } else {
      stop("Latent dimension must be at least 2.")
    }
  }
}
