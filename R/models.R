#' Create a spherical covariance matrix
#'
#' A nimbleFunction that returns a diagonal matrix with all diagonal entries equal to tau_sq. This function is not intended to be called by users of the package directly.
#'
#' @param tau_sq The variance parameter for the spherical covariance matrix.
#' @param p The dimension of the covariance matrix.
#' @return A p x p spherical covariance matrix.
#' @keywords internal
#' @seealso \code{\link[nimble]{nimbleFunction}} for information on nimbleFunctions.
#' @export
makeSphericalSigma <- nimble::nimbleFunction(
  run = function(tau_sq = double(0), p = integer(0)) {
    returnType(double(2))  # returns a p x p matrix
    Sigma <- matrix(0, nrow = p, ncol = p)
    for(i in 1:p){
      Sigma[i,i] <- tau_sq
    }
    return(Sigma)
  }
)

#' Create a diagonal covariance matrix
#'
#' A nimbleFunction that returns a diagonal matrix with diagonal entries equal to tau_sq_vec. This function is not intended to be called by users of the package directly.
#'
#' @param tau_sq_vec A vector of length `p` containing variance parameters for the covariance matrix.
#' @return A p x p spherical covariance matrix.
#' @keywords internal
#' @seealso \code{\link[nimble]{nimbleFunction}} for information on nimbleFunctions.
#' @export
#'
makeDiagonalSigma <- nimble::nimbleFunction(
run = function(tau_sq_vec = double(1)) {
  returnType(double(2))
  p <- length(tau_sq_vec)
  Sigma <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    Sigma[i,i] <- tau_sq_vec[i]
  }
  return(Sigma)
}
)


#' Nimble Code for the Dirichlet Process Clustering with Dissimilarities (DPCD) Models
#'
#' @keywords internal
#' @noRd

unequal_unrestricted_model <- nimble::nimbleCode({
  for (k in 1:(N-1)) {
    beta[k] ~ dbeta(1, alpha_0)
  }
  pi[1:N] <- stick_breaking(beta[1:(N-1)])

  for (k in 1:N) {
    Sigma[1:p,1:p,k] ~ dinvwish(Psi_0[1:p,1:p], nu_0)
    Sigma_mu[1:p,1:p,k] <- (1/lambda) * Sigma[1:p,1:p,k]
    mu[k,1:p] ~ dmnorm(mu_0[1:p], cov = Sigma_mu[1:p,1:p,k])
  }

  sigma_sq ~ dinvgamma(a_0, b_0)

  for (i in 1:n) {
    z[i] ~ dcat(pi[1:N])
    x[i,1:p] ~ dmnorm(mu[z[i], 1:p], cov = Sigma[1:p,1:p,z[i]])
  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      delta[i,j] <- sqrt(sum((x[i,1:p] - x[j,1:p])^2))
      d_obs[i,j] ~ T(dnorm(delta[i,j], sd = sqrt(sigma_sq)), 0, )
    }
  }
})

equal_unrestricted_model <- nimble::nimbleCode({
  for (k in 1:(N-1)) {
    beta[k] ~ dbeta(1, alpha_0)
  }
  pi[1:N] <- stick_breaking(beta[1:(N-1)])

  Sigma[1:p,1:p] ~ dinvwish(Psi_0[1:p,1:p], nu_0)
  Sigma_mu[1:p,1:p] <- (1/lambda) * Sigma[1:p,1:p]

  for (k in 1:N) {
    mu[k,1:p] ~ dmnorm(mu_0[1:p], cov = Sigma_mu[1:p,1:p])
  }

  sigma_sq ~ dinvgamma(a_0, b_0)

  for (i in 1:n) {
    z[i] ~ dcat(pi[1:N])
    x[i,1:p] ~ dmnorm(mu[z[i], 1:p], cov = Sigma[1:p,1:p])
  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      delta[i,j] <- sqrt(sum((x[i,1:p] - x[j,1:p])^2))
      d_obs[i,j] ~ T(dnorm(delta[i,j], sd = sqrt(sigma_sq)), 0, )
    }
  }
})

unequal_diagonal_model <- nimble::nimbleCode({
  for (k in 1:(N-1)) {
    beta[k] ~ dbeta(1, alpha_0)
  }
  pi[1:N] <- stick_breaking(beta[1:(N-1)])

  for (k in 1:N) {
    for (q in 1:p) {
      tau_sq_vec[q, k] ~ dinvgamma(alpha_tau, beta_tau)
    }
    Sigma[1:p,1:p,k] <- makeDiagonalSigma(tau_sq_vec[1:p, k])
    Sigma_mu[1:p,1:p,k] <- (1/lambda) * Sigma[1:p,1:p,k]
    mu[k,1:p] ~ dmnorm(mu_0[1:p], cov = Sigma_mu[1:p,1:p,k])
  }

  sigma_sq ~ dinvgamma(a_0, b_0)

  for (i in 1:n) {
    z[i] ~ dcat(pi[1:N])
    x[i,1:p] ~ dmnorm(mu[z[i], 1:p], cov = Sigma[1:p,1:p,z[i]])
  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      delta[i,j] <- sqrt(sum((x[i,1:p] - x[j,1:p])^2))
      d_obs[i,j] ~ T(dnorm(delta[i,j], sd = sqrt(sigma_sq)), 0, )
    }
  }
})

equal_diagonal_model <- nimble::nimbleCode({
  for (k in 1:(N-1)) {
    beta[k] ~ dbeta(1, alpha_0)
  }
  pi[1:N] <- stick_breaking(beta[1:(N-1)])

  for (q in 1:p) {
    tau_sq_vec[q] ~ dinvgamma(alpha_tau, beta_tau)
  }

  Sigma[1:p,1:p] <- makeDiagonalSigma(tau_sq_vec[1:p])
  Sigma_mu[1:p,1:p] <- (1/lambda) * Sigma[1:p,1:p]
  for (k in 1:N) {
    mu[k,1:p] ~ dmnorm(mu_0[1:p], cov = Sigma_mu[1:p,1:p])
  }

  sigma_sq ~ dinvgamma(a_0, b_0)

  for (i in 1:n) {
    z[i] ~ dcat(pi[1:N])
    x[i,1:p] ~ dmnorm(mu[z[i], 1:p], cov = Sigma[1:p,1:p])
  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      delta[i,j] <- sqrt(sum((x[i,1:p] - x[j,1:p])^2))
      d_obs[i,j] ~ T(dnorm(delta[i,j], sd = sqrt(sigma_sq)), 0, )
    }
  }
})

unequal_spherical_model <- nimble::nimbleCode({
  for (k in 1:(N-1)) {
    beta[k] ~ dbeta(1, alpha_0)
  }
  pi[1:N] <- stick_breaking(beta[1:(N-1)])

  for (k in 1:N) {
    tau_sq[k] ~ dinvgamma(alpha_tau, beta_tau)
    Sigma[1:p,1:p,k] <- makeSphericalSigma(tau_sq[k], p)
    Sigma_mu[1:p,1:p,k] <- (1/lambda) * Sigma[1:p,1:p,k]
    mu[k,1:p] ~ dmnorm(mu_0[1:p], cov = Sigma_mu[1:p,1:p,k])
  }

  sigma_sq ~ dinvgamma(a_0, b_0)

  for (i in 1:n) {
    z[i] ~ dcat(pi[1:N])
    x[i,1:p] ~ dmnorm(mu[z[i], 1:p], cov = Sigma[1:p,1:p,z[i]])
  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      delta[i,j] <- sqrt(sum((x[i,1:p] - x[j,1:p])^2))
      d_obs[i,j] ~ T(dnorm(delta[i,j], sd = sqrt(sigma_sq)), 0, )
    }
  }
})

equal_spherical_model <- nimble::nimbleCode({
  for (k in 1:(N-1)) {
    beta[k] ~ dbeta(1, alpha_0)
  }
  pi[1:N] <- stick_breaking(beta[1:(N-1)])

  tau_sq ~ dinvgamma(alpha_tau, beta_tau)

  Sigma[1:p,1:p] <- makeSphericalSigma(tau_sq, p)
  Sigma_mu[1:p,1:p] <- (1/lambda) * Sigma[1:p,1:p]

  for (k in 1:N) {
    mu[k,1:p] ~ dmnorm(mu_0[1:p], cov = Sigma_mu[1:p,1:p])
  }

  sigma_sq ~ dinvgamma(a_0, b_0)

  for (i in 1:n) {
    z[i] ~ dcat(pi[1:N])
    x[i,1:p] ~ dmnorm(mu[z[i], 1:p], cov = Sigma[1:p,1:p])
  }

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      delta[i,j] <- sqrt(sum((x[i,1:p] - x[j,1:p])^2))
      d_obs[i,j] ~ T(dnorm(delta[i,j], sd = sqrt(sigma_sq)), 0, )
    }
  }
})
