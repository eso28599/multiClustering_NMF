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
    nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat[i][[1]], relErr = "Sil")
    #calculate f difference
    n_views <- length(nmtf_results$Foutput)
    f_score <- c()
    for (j in 1:(n_views - 1)){
        for (k in (j + 1):n_views){
            f_score <- c(f_score,
                 sd(nmtf_results$Foutput[[j]] - nmtf_results$Foutput[[k]]))
        }
    }
    #export each data frame to separate sheets in same Excel file
    file_path <- paste0(data_name, method_vec[i])
    row_nmtf_filename <- paste0(file_path, "_row_clusts.xlsx")
    col_nmtf_filename <- paste0(file_path, "_col_clusts.xlsx")
    error_nmtf_filename <- paste0(file_path, "_errors.csv")
    f_sd_filename <- paste0(file_path, "_f_diff.csv")
    openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
    openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
    suppressMessages(write.csv(mean(f_score), file = f_sd_filename))
    suppressMessages(write.csv(nmtf_results$Error, file = error_nmtf_filename))
}
