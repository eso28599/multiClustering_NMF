#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
row_cl_dims <- list(c(75, 50, 50, 25), c(75, 50, 50, 25))
#100, 50,250 features respectively
col_cl_dims <- list(c(10, 20, 10, 10), c(100, 50, 75, 25))
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

nIter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1,c(2)]<-1


#testing different numbers of clusters
k_vec <- c(2, 3, 4, 5, 6)
test_phi <- c(seq(0,0.02,0.002),0.03,0.04)



#test_phi <- c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 0.65, 0.8, 1.0, 5.0, 10.0, 20.0, 50.0)
