#' Unequal Unrestricted DP Model Code
#'
#' Nimble code for the Unequal Unrestricted model.
#'
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
