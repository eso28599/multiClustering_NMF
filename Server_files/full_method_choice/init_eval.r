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
true_rows <- import_matrix(paste0(data_name, "/true_rows.xlsx"))
true_cols <- import_matrix(paste0(data_name, "/true_cols.xlsx"))

data_views <- import_matrix(paste0(data_name, "/data.xlsx"), col_n = TRUE)
end_vec <- c("Iter", "Conv")
choice_vec <- c("Error", "RelErr", "Sil")

for (end in end_vec){
  for (choice in choice_vec){
    file_path <- paste0(data_name, paste0("/", paste0(end, choice)))
    row_filename <- paste0(file_path, "_nmtf_row_clusts.xlsx")
    col_filename <- paste0(file_path, "_nmtf_col_clusts.xlsx")
    error_filename <- paste0(file_path, "_nmtf_errors.csv")
    #results for nmtf
    results_nmtf <- evaluate_simulation(import_matrix(row_filename),
                                              import_matrix(col_filename),
                                              true_rows, true_cols, data_views,
                                              read.csv(error_filename, row.names = 1))
    #export each data frame to separate sheets in same Excel file
    openxlsx::write.xlsx(results_nmtf, file = paste0(file_path, "_results_nmtf.xlsx"))
  }
}
