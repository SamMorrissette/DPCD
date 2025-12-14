set.seed(123)

# Generate samples from a two-component mixture model (n = 10)
x <- rbind(matrix(rnorm(20), ncol = 2), matrix(rnorm(20, mean = 4), ncol = 2))
dis_mat <- dist(x)

# Fit the Equal Spherical DPCD model
fit <- run_dpcd(dis_mat, model = "ES", p = 2, niter = 5000, nburn = 1000)
mcmc_example <- fit$samples

usethis::use_data(mcmc_example, overwrite = TRUE)
