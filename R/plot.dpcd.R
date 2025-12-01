plot.dpcd.2D <- function(mcmc_samples, align_matrix) {
  if (!inherits(mcmc_samples, "mcmc.list") & !inherits(mcmc_samples, "mcmc")) {
    stop("Input must be an object of class 'mcmc.list' or 'mcmc'.")
  }
  if (inherits(mcmc_samples, "mcmc.list")) {
    full_samples <- do.call(rbind, mcmc_samples)
  } else {
    full_samples <- mcmc_samples
  }

  x_samples <- full_samples[, startsWith(colnames(full_samples), "x")]
  z_samples <- full_samples[, startsWith(colnames(full_samples), "z")]
  pi_samples <- full_samples[, startsWith(colnames(full_samples), "pi")]
  mu_samples <- full_samples[, startsWith(colnames(full_samples), "mu")]
  Sigma_samples <- full_samples[, startsWith(colnames(full_samples), "Sigma")]

  iters <- nrow(full_samples)
  x_transformed <- matrix(NA, nrow = iters, ncol = ncol(x_samples))
  for (i in 1:iters) {
    x_mat <- matrix(x_samples[i,], ncol = 2, byrow = FALSE)
    x_transformed[i,] <- c(procrustes(align_matrix, x_mat))
  }
  final_mat <- matrix(colMeans(x_transformed), ncol = 2, byrow = FALSE)


  psm <- mcclust::comp.psm(z_samples)
  labels <- mcclust::maxpear(psm)
  plot(final_mat, col = labels$cl, pch = 19)
}
