source("extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
#row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
row_cl_dims <- list(rep(50, 3), rep(200, 3), rep(300, 3), rep(500, 3), rep(1000, 3))

#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 50 individuals 75:75:50 -> 1.5:1.5: 17:17:12
row_start_50 <- list(c(1, 19, 39), c(1, 19, 39), c(1, 19, 39))
row_end_50 <- list(c(18, 38, 50), c(18, 38, 50), c(18, 38, 50))
# 200 individuals
row_start_200 <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end_200 <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
# 300 individuals + 35, 35, 30
row_start_300 <- list(c(1, 111, 221), c(1, 111, 221), c(1, 111, 221))
row_end_300 <- list(c(110, 220, 300), c(110, 220, 300), c(110, 220, 300))
# 500 individuals + 35, 35, 30
row_start_500 <- list(c(1, 186, 371), c(1, 186, 371), c(1, 186, 371))
row_end_500 <- list(c(185, 370, 500), c(185, 370, 500), c(185, 370, 500))
# 1000 individuals + 15, 35, 30
row_start_1000 <- list(c(1, 371, 741), c(1, 371, 741), c(1, 371, 741))
row_end_1000 <- list(c(370, 740, 1000), c(370, 740, 1000), c(370, 740, 1000))
row_start <- list(row_start_50, row_start_200, row_start_300,
                 row_start_500, row_start_1000)
row_end <- list(row_end_50, row_end_200, row_end_300,
                     row_end_500, row_end_1000)
col_start <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151))
col_end <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250))

niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
val <- 200
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
phi_mat <- val*phi
#phi_mat <- 0.01*phi
method_vec <- c("/res_nmtf_50", "/res_nmtf_200", "/res_nmtf_300", "/res_nmtf_500", "/res_nmtf_1000")

numbers <- c(50, 200, 300, 500, 1000)
method_vec_sgl <- paste0("/nmtf_", numbers)
method_vec_gfa <- paste0("/gfa_", numbers)
method_vec_issvd <- paste0("/issvd_", numbers)
factor <- "indiv"
kept_factor <- paste0(length(col_start), " views, ", 
            length(col_start[1]), " biclusters, ")
factor_name <- "Individuals"
x_title <- "Number of individuals"
plot_title <- "The effect of increasing the number of inidividuals"
file_names <- c(paste0("res_nmtf_", numbers),
    paste0("gfa_", numbers), 
    paste0("issvd_", numbers), paste0("nmtf_", numbers))
method_vec_res <- c(rep("ResNMTF", 5), rep("GFA", 5),
 rep("iSSVD", 5), rep("NMTF", 5))
k_vec <- rep(3, length=length(method_vec_res))
factor_vec <- rep(numbers, 4)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- numbers
