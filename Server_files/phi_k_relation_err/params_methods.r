args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1]) #simulation folder
batch_folder = as.character(args[2]) #repeat
library("curl")
library("readxl")
library("rio")
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("resNMTF_funcs.r")
#import as list of matrices instead of as a list
import_matrix <- function(filename){
    return(lapply(import_list(filename), function(x) as.matrix(x)))
}
#read in data
data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
data_views <- import_matrix(paste0(data_name, "/data.xlsx"))
#now apply for each value of phi we want to consider
for (k in k_vec){
  no_row_cl <- rep(k, n_views)
  no_col_cl <- rep(k, n_views)
  for (phi_val in test_phi){
    phi_mat <- phi_val*phi
    nmtf_results <- restMultiNMTF_run(Xinput = data_views, KK = no_row_cl,
     LL = no_col_cl, phi = phi_mat, nIter = nIter)
    #export each data frame to separate sheets in same Excel file
    file_path <- paste0(data_name, paste0("/", paste0(k, phi_val)))
    row_nmtf_filename <- paste0(file_path, "_nmtf_row_clusts.xlsx")
    col_nmtf_filename <- paste0(file_path, "_nmtf_col_clusts.xlsx")
    error_nmtf_filename <- paste0(file_path, "_nmtf_errors.csv")
    openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
    openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
    write.csv(nmtf_results$Error, file = error_nmtf_filename)
  }
}
