#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
#row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
row_cl_dims <- list(rep(50, 3), rep(200, 3), rep(300, 3))

#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 50 individuals 75:75:50 -> 1.5:1.5: 17:17:12
row_start_50 <- list(c(1, 19, 39), c(1, 19, 39), c(1, 19, 39))
row_end_50 <- list(c(18, 38, 50), c(18, 38, 50), c(18, 38, 50))
# 200 individuals
row_start_200 <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end_200 <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
# 300 individuals + 35, 35, 30
row_start_300 <- list(c(36, 146, 251), c(36, 146, 251), c(36, 146, 251))
row_end_300 <- list(c(110, 220, 300), c(110, 220, 300), c(110, 220, 300))
col_start_200 <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151))
col_end_200 <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250))
row_start <- list(row_start_50, row_start_200, row_start_300)
row_end <- list(row_end_50, row_end_200, row_end_300)

col_start <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151))
col_end <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250))

niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
phi_mat <- phi
#phi_mat <- 0.01*phi
method_vec <- c("/res_nmtf_50", "/res_nmtf_200", "/res_nmtf_300")
