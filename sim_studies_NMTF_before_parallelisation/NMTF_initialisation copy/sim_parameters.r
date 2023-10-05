#file is specific to simulation
#simulation parameters
n_repeats <- 3
noise_level <- 3
# exhausitive, exclusive clusters
#each dataset with 200 samples, same clusters across views
row_cl_dims <- list(c(75, 75, 50),c(75, 75, 50),c(75, 75, 50))
#100, 50,250 features respectively
col_cl_dims <- list(c(30,30,40),c(10, 20, 20),c(100,50,100))

# where are datasets stored - from current location on server - should just be /sim_name_data
path_to_data_folder <- "sim_studies_NMTF/sim1_data"
#number of datasets
# can we just count number of files or does this have to be inputed
n_views <- 3
no_row_cl <- c(3,3,3)
no_col_cl <- c(3,3,3)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1,2:3]<-0.5
nIter <- 1000
