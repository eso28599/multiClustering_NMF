#file is specific to simulation
#simulation parameters

#each dataset with 200 samples, same clusters across views
#row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
row_cl_dims <- rep(200, 3)

#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 3 biclusters
row_start <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
col_start <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151))
col_end <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250))

niter <- 1000
noise_level <- 1:10
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
method_vec <- paste0("/res_nmtf_", noise_level)
