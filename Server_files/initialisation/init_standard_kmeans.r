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
#apply method
nmtf_results <- restMultiNMTF_run(Xinput = data_views, KK = no_row_cl,
     LL = no_col_cl, phi = phi, nIter = nIter)

#export each data frame to separate sheets in same Excel file
row_nmtf_filename <- paste0(data_name, "/nmtf_row_clusts.xlsx")
col_nmtf_filename <- paste0(data_name, "/nmtf_col_clusts.xlsx")
error_nmtf_filename <- paste0(data_name, "/nmtf_errors.xlsx")
openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
openxlsx::write.xlsx(nmtf_results$Error, file = error_nmtf_filename)


#k- means initialisation
source(paste0(path_to_sim_folder, "/resNMTF_k_means.r"))
#apply method
nmtf_results <- restMultiNMTF_run(Xinput = data_views, KK = no_row_cl,
     LL = no_col_cl, phi = phi, nIter = nIter)

#export each data frame to separate sheets in same Excel file
row_nmtf_filename <- paste0(data_name, "/nmtf_k_means_row_clusts.xlsx")
col_nmtf_filename <- paste0(data_name, "/nmtf_k_means_col_clusts.xlsx")
error_nmtf_filename <- paste0(data_name, "/nmtf_k_means_errors.xlsx")
openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
openxlsx::write.xlsx(nmtf_results$Error, file = error_nmtf_filename)