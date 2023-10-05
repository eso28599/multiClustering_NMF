#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
row_cl_dims <- list(c(75, 50, 50, 25), c(75, 50, 50, 25))
#100, 50,250 features respectively
col_cl_dims <- list(c(10, 20, 10, 10), c(100, 50, 75, 25))
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, 2] <- 1
phi_mat <- 0.01*phi
