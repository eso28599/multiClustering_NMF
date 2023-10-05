#file is specific to simulation
#simulation parameters
noise_level <- 3
#each dataset with 200 samples, same clusters across views
row_cl_dims <- list(c(75, 75, 50),c(75, 75, 50),c(75, 75, 50))
#100, 50,250 features respectively
col_cl_dims <- list(c(30,30,40),c(10, 20, 20),c(100,50,100))
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

nIter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1,c(2,3)]<-1
phi[2,c(3)]<-1

test_phi <- seq(50,100,10)
