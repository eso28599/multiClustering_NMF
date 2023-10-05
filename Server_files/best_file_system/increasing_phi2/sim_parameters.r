#file is specific to simulation
#simulation parameters

#each dataset with 200 samples, same clusters across views
#row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
row_cl_dims <- rep(200, 4)

#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250, 150)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 2 biclusters
row_start <- list(c(1, 76), c(1, 76), c(1, 76), c(1, 76))
row_end <- list(c(75, 200), c(75, 200),
            c(75, 200), c(75, 200))
col_start <- list(c(1, 31), c(1, 11), c(1, 101), c(1, 51))
col_end <- list(c(30, 100), c(10, 50),
                 c(100, 250), c(50, 150))

niter <- 1000
noise_level <- 5
# exhausitive, exclusive cluster
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi_vals <- c(1, seq(200, 3000, 400), seq(6000,50000, 2000))
phi[1, c(2, 3, 4)] <- 1
phi[2, c(3, 4)] <- 1
phi[3, 4] <- 1
#phi_mat <-  vector("list", length = length(phi_vals))
phi_mat <- lapply(phi_vals, function(x) x * phi)
#phi_mat <- 0.01*phi
method_vec <- paste0("/res_nmtf_", phi_vals)
