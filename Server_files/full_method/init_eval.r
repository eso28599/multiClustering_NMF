args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
batch_folder = as.character(args[2])

# simulation evaluation file
source(paste0(path_to_sim_folder, "/sim_parameters.r")) #parameters for simulation
source("evaluation_funcs.r")
library("rio")
#import as list of matrices instead of as a list
import_matrix <- function(filename, col_n = FALSE){
    return(lapply(import_list(filename, col_names=col_n), function(x) as.matrix(x)))
}
#read in results
data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
#export each data frame to separate sheets in same Excel file
true_row_filename <- paste0(data_name, "/true_rows.xlsx")
true_col_filename <- paste0(data_name, "/true_cols.xlsx")
row_nmtf_filename <- paste0(data_name, "/nmtf_row_clusts.xlsx")
col_nmtf_filename <- paste0(data_name, "/nmtf_col_clusts.xlsx")
row_k_means_nmtf_filename <- paste0(data_name, "/nmtf_k_means_row_clusts.xlsx")
col_k_means_nmtf_filename <- paste0(data_name, "/nmtf_k_means_col_clusts.xlsx")
data_views <- import_matrix(paste0(data_name, "/data.xlsx"), col_n = TRUE)
error_filename <- paste0(data_name, "/nmtf_errors.csv")
error_k_means_filename <- paste0(data_name, "/nmtf_k_means_errors.csv")

#results for nmtf
results_nmtf <- evaluate_simulation(import_matrix(row_nmtf_filename),
                                            import_matrix(col_nmtf_filename),
                                        import_matrix(true_row_filename),
                                        import_matrix(true_col_filename),
                                         data_views,
                                        read.csv(error_filename, row.names = 1))
results_k_means_nmtf <- evaluate_simulation(import_matrix(row_k_means_nmtf_filename),
                                        import_matrix(col_k_means_nmtf_filename),
                                        import_matrix(true_row_filename),
                                        import_matrix(true_col_filename),
                                         data_views,
                                        read.csv(error_k_means_filename, row.names = 1))


#save results as worksheets
#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(results_nmtf, file = paste0(data_name, "/results_nmtf.xlsx"))
openxlsx::write.xlsx(results_k_means_nmtf, file = paste0(data_name, "/results_k_means_nmtf.xlsx"))
