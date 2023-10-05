args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1]) #simulation folder
batch_folder = as.character(args[2]) #repeat
suppressPackageStartupMessages(library("curl"))
library("readxl")
suppressPackageStartupMessages(library("rio"))
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
source("resNMTF_funcs.r")
#import as list of matrices instead of as a list
import_matrix <- function(filename, col_n = FALSE){
    return(lapply(import_list(filename, col_names=col_n), function(x) as.matrix(x)))
}
#read in data
data_name <- paste0(paste0(path_to_sim_folder, "/data/"), batch_folder)
data_views <- import_matrix(paste0(data_name, "/data.xlsx"), col_n = TRUE)

for (i in 1:length(method_vec)){
    #nmtf
    set.seed(as.numeric(batch_folder) + i)
    data_sub_index <- sample(length(data_views), i+1)
    data_subset <- data_views[data_sub_index]
    nmtf_results <- restMultiNMTF_run(Xinput = data_subset, phi = phi_mat[i][[1]], relErr = "Sil")
    #export each data frame to separate sheets in same Excel file
    file_path <- paste0(data_name, method_vec[i])
    row_nmtf_filename <- paste0(file_path, "_row_clusts.xlsx")
    col_nmtf_filename <- paste0(file_path, "_col_clusts.xlsx")
    error_nmtf_filename <- paste0(file_path, "_errors.csv")
    openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
    openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
    suppressMessages(write.csv(nmtf_results$Error, file = error_nmtf_filename))
}
