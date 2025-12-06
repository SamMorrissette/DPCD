#' Procrustes Transformation
#'
#' Aligns a given object configuration to a target object configuration using a Procrustes transformation.
#'
#' @param X The target configuration.
#' @param Y The configuration to be aligned to X.
#'
#' @details This function performs a Procrustes transformation to align a given configuration, `Y`, to the target configuration, `X`, using a combination of translation and rotation. The transformation aims to minimize the sum of squared differences between the two configurations.
#'
#' `X` and `Y` should be numeric matrices of the same dimension.
#'
#' @return The transformed version of Y aligned to X.
#'
#' @examples
#' X <- matrix(rnorm(20), ncol = 2)
#' rotation_matrix <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)), ncol = 2)
#' Y <- X %*% rotation_matrix + 2
#' Y_transformed <- procrustes(X, Y)
#' @export

procrustes <- function(X, Y) {
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)

  C <- t(X_centered) %*% Y_centered

  svd_result <- svd(C)

  R <- svd_result$v %*% t(svd_result$u)

  Y_rotated <- Y_centered %*% R

  t_vec = colMeans(X) - colMeans(Y_rotated)
  Y_transformed <- sweep(Y_rotated, 2, t_vec, "+")

  return(Y_transformed)
}
