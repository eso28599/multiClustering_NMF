#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
#row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
row_cl_dims <- rep(200, 4)
row_start <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250, 150)
col_start <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151), c(1, 41, 91))
col_end <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250), c(40, 90, 150))
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3, 4)] <- 1
phi[2, c(3, 4)] <- 1
phi[3, c(4)] <- 1
phi_mat_2 <- phi[1:2, 1:2]
phi_mat_3 <- phi[1:3, 1:3]
phi_mat_4 <- phi
phi_mat <- list(phi_mat_2, phi_mat_3, phi_mat_4)
#phi_mat <- 0.01*phi
method_vec <- c("/res_nmtf_2", "/res_nmtf_3", "/res_nmtf_4")
