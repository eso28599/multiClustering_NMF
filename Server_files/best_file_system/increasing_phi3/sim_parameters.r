#file is specific to simulation
#simulation parameters

#each dataset with 200 samples, same clusters across views
row_cl_dims <- rep(200, 2)

#100, 50,250 features respectively
col_cl_dims <- c(100, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

#2 views, 4 biclusters
row_start <- list(c(1, 76, 111, 151), c(1, 76, 111, 151))
row_end <- list(c(75, 110, 150, 200), c(75, 110, 150, 200))
col_start <- list(c(1, 31, 46, 61),  c(1, 101, 126, 151))
col_end <- list(c(30, 45, 60, 100), c(100, 125, 150, 250))

niter <- 1000
noise_level <- 5
# exhausitive, exclusive cluster
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi_vals <- c(1, seq(200, 3000, 400), seq(6000,50000, 2000))
phi[1, 2] <- 1
phi_mat <- lapply(phi_vals, function(x) x * phi)
method_vec <- paste0("/res_nmtf_", phi_vals)
