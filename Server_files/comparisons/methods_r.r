args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1]) #simulation folder
batch_folder = as.character(args[2]) #repeat
suppressPackageStartupMessages(library("curl"))
library("readxl")
library("rio")
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("resNMTF_funcs.r")
#import as list of matrices instead of as a list
import_matrix <- function(filename, col_n = FALSE){
    return(lapply(import_list(filename, col_names=col_n), function(x) as.matrix(x)))
}
#read in data
data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
data_views <- import_matrix(paste0(data_name, "/data.xlsx"), col_n = TRUE)

#nmtf
nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat, relErr = "Sil")
#export each data frame to separate sheets in same Excel file
file_path <- paste0(data_name, "/res")
row_nmtf_filename <- paste0(file_path, "_nmtf_row_clusts.xlsx")
col_nmtf_filename <- paste0(file_path, "_nmtf_col_clusts.xlsx")
error_nmtf_filename <- paste0(file_path, "_nmtf_errors.csv")
openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
suppressMessages(write.csv(nmtf_results$Error, file = error_nmtf_filename))


#nmtf_results <- restMultiNMTF_run_mix(Xinput = data_views, phi = phi_mat, relErr = "RelErr")
#export each data frame to separate sheets in same Excel file
#file_path <- paste0(data_name, "/mix")

#row_nmtf_filename <- paste0(file_path, "_nmtf_row_clusts.xlsx")
#col_nmtf_filename <- paste0(file_path, "_nmtf_col_clusts.xlsx")
#error_nmtf_filename <- paste0(file_path, "_nmtf_errors.csv")
#openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
#openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
#suppressMessages(write.csv(nmtf_results$Error, file = error_nmtf_filename))



#row_mvbc_filename <- paste0(data_name, "/mvbc_row_clusts.xlsx")
#col_mvbc_filename <- paste0(data_name, "/mvbc_col_clusts.xlsx")
#openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_mvbc_filename)
#openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_mvbc_filename)