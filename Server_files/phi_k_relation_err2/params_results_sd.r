args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
folder_list = args[2]
source("summary_funcs.r")
source(paste0(path_to_sim_folder, "/sim_parameters.r")) #parameters for simulation
library(dplyr)
library(kableExtra)

#file to extract results
#read in folder_list
folder_names <- read.delim(folder_list, colClasses = "character", header=FALSE)[,1]

n_repeats <- length(folder_names)
# n_repeats x 3 matrices
#measures
#measures <- c("Accuracy", "ARI", "NMI", "Sim", "Dist", "Err")
measures <- c("ARI", "NMI", "Sil", "Error")
#averages <- c("Row", "Column", "Overall")
averages <- c("Overall")

#initialise storage of results
#average col and row values, as well as total average
#now want to get averages and standard deviations across repetitions
#store results in a dataframe

results <- data.frame(Method=c(), Cluster=c(), Phi=c(), Mean=c(), Sd=c())
for (k in k_vec){
  no_row_cl <- rep(k, n_views)
  no_col_cl <- rep(k, n_views)
  avgs_nmtf <- replicate(length(test_phi), matrix(nrow = n_repeats, ncol = length(averages)), simplify=FALSE)
  for (phi_val in test_phi){
    for (i in 1:length(folder_names)){
        data_name <- paste0(paste0(path_to_sim_folder, "/data/"), folder_names[[i]])
        avgs_nmtf <- get_averages(avgs_nmtf, data_name, paste0(paste0(k, phi_val),"_results_nmtf"), measures, i)
    }
    for (i in 1:length(measures)){
        results <- rbind(results,
                          avg_over_repeats(avgs_nmtf,i,measures[i],k,phi_val))
        }
    }
    #store results - one workbook for each test - each sheet a different measure
    #save_file <- paste0(paste0(path_to_sim_folder,paste0("/",paste0(k, phi_val))), "_avgs_nmtf.xlsx")
    #openxlsx::write.xlsx(avgs_nmtf, file = save_file)
  }
#list of files

#change mean and sd columns to numeric 
colnames(results) <- c("k", "Phi", "Measure", "Mean", "S.d.")
results <- mutate_at(results,c("Mean","S.d."), as.numeric)
write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))

#make and saves plots
path_to_save <- paste0(paste0(path_to_sim_folder,"/"), paste0("Overall_diff_k_and_phi.pdf"))
make_plot_phi(subset(results,measures=="Error"), "Overall performance with different values of k",
                paste0(n_repeats, " repetitions"),
                path_to_save)

path_to_save <- paste0(path_to_sim_folder, "/diff_params_table_")
#produce tables
#all_table(results,test_phi,path_to_save, feat="k",group=FALSE)
all_table2(results,k_vec,test_phi, paste0(path_to_save,2), feat="Phi",group=FALSE)