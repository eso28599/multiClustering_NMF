source("extra_funcs.r")
#file is specific to simulation
#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
#row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
row_cl_dims <- rep(200, 3)

#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 2 biclusters
row_start_2 <- list(c(1, 76), c(1, 76), c(1, 76))
row_end_2 <- list(c(75, 200), c(75, 200), c(75, 200))
col_start_2 <- list(c(1, 31), c(1, 11), c(1, 101))
col_end_2 <- list(c(30, 100), c(10, 50), c(100, 250))
# 3 biclusters
row_start_3 <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end_3 <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
col_start_3 <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151))
col_end_3 <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250))
# 4 biclusters
row_start_4 <- list(c(1, 76, 111, 151), c(1, 76, 111, 151), c(1, 76, 111, 151))
row_end_4 <- list(c(75, 110, 150, 200), c(75, 110, 150, 200),
                 c(75, 110, 150, 200))
col_start_4 <- list(c(1, 31, 46, 61), c(1, 11, 21, 31), c(1, 101, 126, 151))
col_end_4 <- list(c(30, 45, 60, 100), c(10, 20, 30, 50), c(100, 125, 150, 250))
# 5 biclusters
row_start_5 <- list(c(1, 51, 76, 111, 151), c(1, 51, 76, 111, 151),
                     c(1, 51, 76, 111, 151))
row_end_5 <- list(c(50, 75, 110, 150, 200), c(50, 75, 110, 150, 200),
                 c(50, 75, 110, 150, 200))
col_start_5 <- list(c(1, 21, 31, 46, 61), c(1, 6, 11, 21, 31),
                         c(1, 51, 101, 126, 151))
col_end_5 <- list(c(20, 30, 45, 60, 100), c(5, 10, 20, 30, 50),
                     c(50, 100, 125, 150, 250))
# 6 biclusters
row_start_6 <- list(c(1, 51, 76, 111, 151, 176), c(1, 51, 76, 111, 151, 176),
                     c(1, 51, 76, 111, 151, 176), c(1, 51, 76, 111, 151, 176))
row_end_6 <- list(c(50, 75, 110, 150, 175, 200), c(50, 75, 110, 150, 175, 200),
                 c(50, 75, 110, 150, 175, 200), c(50, 75, 110, 150, 175, 200))
col_start_6 <- list(c(1, 21, 31, 46, 61, 81), c(1, 6, 11, 21, 31, 41),
                         c(1, 51, 101, 126, 151, 201), c(1, 26, 51, 76, 101, 126))
col_end_6 <- list(c(20, 30, 45, 60, 80, 100), c(5, 10, 20, 30, 40, 50),
                     c(50, 100, 125, 150, 200, 250),
                      c(25, 50, 75, 100, 125, 150))
row_start <- list(row_start_2, row_start_3,
             row_start_4, row_start_5, row_start_6)
row_end <- list(row_end_2, row_end_3, row_end_4, row_end_5, row_end_6)
col_start <- list(col_start_2, col_start_3,
             col_start_4, col_start_5, col_start_6)
col_end <- list(col_end_2, col_end_3, col_end_4, col_end_5, col_end_6)

niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3)] <- 1
phi[2, c(3)] <- 1
val <- 200
phi_mat <- val * phi
bicl_numbers <- 2:6
method_vec <- c("/res_nmtf_2", "/res_nmtf_3",
 "/res_nmtf_4", "/res_nmtf_5", "/res_nmtf_6")
method_vec_sgl <- paste0("/nmtf_", bicl_numbers)
method_vec_gfa <- paste0("/gfa_", bicl_numbers)
method_vec_issvd <- paste0("/issvd_", bicl_numbers)
factor <- "bicl"
kept_factor <- paste0(length(row_cl_dims), " views, ")
factor_name <- "Biclusters"
x_title <- "Number of biclusters"
plot_title <- "The effect of increasing the number of biclusters on performance"
file_names <- c(paste0("res_nmtf_", bicl_numbers),
    paste0("gfa_", bicl_numbers), 
    paste0("issvd_", bicl_numbers), paste0("nmtf_", bicl_numbers))
method_vec_res <- c(rep("ResNMTF", 5), rep("GFA", 5),
 rep("iSSVD", 5), rep("NMTF", 5))
factor_vec <- rep(2:6, 4)
k_vec <- rep(2:6, 4)
phi_constant <- TRUE
n_col_plots <- 2
n_row_plots <- 2
col_names_tables <- 2:6