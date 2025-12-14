set.seed(123)

# Generate samples from a two-component mixture model (n = 50)
x <- rbind(matrix(rnorm(20), ncol = 2), matrix(rnorm(30, mean = 4), ncol = 2))
dis_mat <- dist(x)

# Fit the Equal Spherical DPCD model
fit <- run_dpcd(dis_mat, model = "ES", p = 2, niter = 50000, nburn = 10000)

usethis::use_data(mcmc_example, overwrite = TRUE)
