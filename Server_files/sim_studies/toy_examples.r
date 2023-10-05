source("sim_studies/data_generation_toy.r")
source("sim_studies/evaluation_funcs.r")
source("sim_studies/resNMTF_funcs.r")
#toy examples
#a -  same biclusters in terms of columns across two of the views
##view 2 column dimensions becomes 100: (30, 30, 40)#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
row_dims <-  rep(200, 3)
#100, 100,250 features respectively
col_dims <- c(100, 100, 250)

row_start <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
col_start <- list(c(1, 31, 61), c(1, 31, 61), c(1, 101, 151))
col_end <- list(c(30, 60, 100),c(30, 60, 100), c(100, 150, 250))
niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_dims)
no_row_cl <- sapply(row_dims, length)
no_col_cl <- sapply(col_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
res_val <- 50000
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
phi_mat <- res_val * phi
#g restrictions
psi <- matrix(0, nrow = n_views, ncol = n_views)
psi[1, c(2)] <- 1
psi_mat <- res_val * psi


#generate data
data_all <- multi_view(row_dims,
        col_dims, row_start, row_end, col_start, col_end, 
        noise=noise_level,row_same_shuffle=TRUE, col_same_shuffle=c(T, T, F), seed=10)
source("sim_studies/increasing_bicl/resNMTF_funcs.r")
nmtf_results <- restMultiNMTF_run(Xinput = data_all$data_views, 
            phi = phi_mat, psi = psi_mat, relErr = "Sil")
#evaluate results
source("sim_studies/evaluation_funcs.r")
nmtf_scores <- evaluate_simulation_comp(nmtf_results$row_clusters , nmtf_results$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)
