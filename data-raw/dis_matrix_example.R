set.seed(123)

# Generate samples from a two-component mixture model (n = 20)
x <- rbind(matrix(rnorm(20), ncol = 2), matrix(rnorm(20, mean = 4), ncol = 2))
dis_mat_example <- dist(x)

usethis::use_data(dis_mat_example, overwrite = TRUE)
