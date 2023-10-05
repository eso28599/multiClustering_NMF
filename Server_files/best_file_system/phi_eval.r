args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
batch_folder = as.character(args[2])

# simulation evaluation file
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("evaluation_funcs.r")
suppressPackageStartupMessages(library("rio"))
#import as list of matrices instead of as a list
import_matrix <- function(filename, col_n = TRUE){
    return(suppressMessages(lapply(import_list(filename,
            col_names = col_n), function(x) as.matrix(x))))
}
#read in results
data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
#export each data frame to separate sheets in same Excel file
true_rows <- import_matrix(paste0(data_name, "/true_rows.xlsx"))
true_cols <- import_matrix(paste0(data_name, "/true_cols.xlsx"))
data_views <- import_matrix(paste0(data_name, "/data.xlsx"), col_n = TRUE)
tol <- 1e-04
#analyse each method
for (method in method_vec){
    file_path <- paste0(data_name, method)
    row_filename <- paste0(file_path, "_row_clusts.xlsx")
    col_filename <- paste0(file_path, "_col_clusts.xlsx")
    f_sd_filename <- paste0(file_path, "_f_diff.csv")
    #results for nmtf
    results <- evaluate_simulation_comp(import_matrix(row_filename),
                                              import_matrix(col_filename),
                                              true_rows,
                                               true_cols,
                                                data_views)
    f_sd_mat <- matrix(as.numeric(read.csv(f_sd_filename)$x < tol),
                 2, dim(results$CSR)[2])
    rownames(f_sd_mat) <- c("Row-clustering", "Column-clustering")
    results[["Restrictions"]] <- f_sd_mat
    #export each data frame to separate sheets in same Excel file
    path_to_save <- paste0(file_path, "_results.xlsx")
    openxlsx::write.xlsx(results, file = path_to_save)
}
# issvd
# results_rep <- vector(mode = "list", length = 5)
# sil_list <- c()
# for (i in 1:5){
#     file_path <- paste0(data_name, paste0("/issvd", (i - 1)))
#     row_filename <- paste0(file_path, "_row_clusts.xlsx")
#     col_filename <- paste0(file_path, "_col_clusts.xlsx")
#     #results for nmtf
#     results_rep[[i]] <- evaluate_simulation_comp(import_matrix(row_filename),
#                                               import_matrix(col_filename),
#                                               true_rows, true_cols, data_views)
#     sil_list <- c(sil_list, mean(results_rep[[i]]$BiS))
# }
# path_to_save <- paste0(data_name, "/issvd_results.xlsx")
# openxlsx::write.xlsx(results_rep[[which.max(sil_list)]], file = path_to_save)
