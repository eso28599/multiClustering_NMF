# this file should be common to all simulations in set up
# but changed to carry out specific simulation
#biclustering approaches in R
source("resNMTF.r")
source("sim_parameters.r")
library("readxl")
#change to just sim_parameters
# write list of names of data views
dataset_names <- sapply(1:n_repeats, function(x) {
    paste0(paste0(paste0(path_to_data_folder, "/data/repeat"), x), ".xlsx")
})
# initialise storage of results
row_cl_returned <- vector("list", length = n_repeats)
col_cl_returned <- vector("list", length = n_repeats)
errors_nmtf <- vector("list", length = n_repeats)
# read in data files and apply relevant methods
for (i in 1:n_repeats){
    data_views <- import_list(dataset_names[i]) #read in data
    #apply restrictive NMTF with SVD based initialisation
    nmtf_results <- restMultiNMTF_run(Xinput = data_views, KK = no_row_cl,
     LL = no_col_cl, phi = phi, nIter = nIter)
    #save results as worksheets
    #export each data frame to separate sheets in same Excel file
    row_nmtf_filename <- paste0(paste0(paste0(path_to_data_folder,'/nmtf_row_clusts/repeat'),i),'.xlsx')
    col_nmtf_filename <- paste0(paste0(paste0(path_to_data_folder,'/nmtf_col_clusts/repeat'),i),'.xlsx')
    error_nmtf_filename <- paste0(paste0(paste0(path_to_data_folder,'/nmtf_errors/repeat'),i),'.xlsx')
    openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
    openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
    openxlsx::write.xlsx(nmtf_results$Error, file = error_nmtf_filename)

    #now apply with k means based initialisation

    nmtf_results <- restMultiNMTF_run(Xinput = data_views, KK = no_row_cl,
     LL = no_col_cl, phi = phi, nIter = nIter)
    #save results as worksheets
    #export each data frame to separate sheets in same Excel file
    row_nmtf_filename <- paste0(paste0(paste0(path_to_data_folder,'/nmtf_row_clusts/repeat'),i),'.xlsx')
    col_nmtf_filename <- paste0(paste0(paste0(path_to_data_folder,'/nmtf_col_clusts/repeat'),i),'.xlsx')
    error_nmtf_filename <- paste0(paste0(paste0(path_to_data_folder,'/nmtf_errors/repeat'),i),'.xlsx')
    openxlsx::write.xlsx(nmtf_results$row_clusters, file = row_nmtf_filename)
    openxlsx::write.xlsx(nmtf_results$col_clusters, file = col_nmtf_filename)
    openxlsx::write.xlsx(nmtf_results$Error, file = error_nmtf_filename)
}
