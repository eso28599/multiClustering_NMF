source("extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
row_cl_dims <- rep(100, 2)

#100, 50,250 features respectively
col_cl_dims <- c(1000, 1000)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2)] <- 1
val <- 200
phi_mat <- val * phi

numbers <- c(1, 2, 5, 10)
method_vec <- paste0("/res_nmtf_", numbers)
method_vec_sgl <- paste0("/nmtf_", numbers)
method_vec_gfa <- paste0("/gfa_", numbers)
method_vec_issvd <- paste0("/issvd_", numbers)
factor <- "ssvd"
kept_factor <- paste0(length(col_cl_dims), " views, ", 
           4, " biclusters, ")
factor_name <- "Scalar"
x_title <- "Scalar"
plot_title <- "The effect of increasing the scale factor"
file_names <- c(paste0("res_nmtf_", numbers),
    paste0("gfa_", numbers), 
    paste0("issvd_", numbers), paste0("nmtf_", numbers))
method_vec_res <- c(rep("ResNMTF", 4), rep("GFA", 4),
 rep("iSSVD", 4), rep("NMTF", 4))
factor_vec <- rep(numbers, 4)
k_vec <- rep(4, 4*4)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- numbers
