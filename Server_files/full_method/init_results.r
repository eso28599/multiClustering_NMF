args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
folder_list = args[2]
source("summary_funcs.r")
library(dplyr)
library(kableExtra)

#file to extract results
#read in folder_list
folder_names <- read.delim(folder_list, colClasses = "character", header=FALSE)[,1]

n_repeats <- length(folder_names)
# n_repeats x 3 matrices
#measures
measures <- c("ARI", "NMI", "Sil", "Error")
averages <- c("Overall")

#initialise storage of results
avgs_nmtf <- replicate(length(measures), matrix(nrow = n_repeats, ncol = length(averages)), simplify=FALSE)
avgs_nmtf_k_means <- replicate(length(measures), matrix(nrow = n_repeats, ncol = length(averages)), simplify=FALSE)
for (i in 1:length(folder_names)){
    data_name <- paste0(paste0(path_to_sim_folder, "/data/"), folder_names[[i]])
    avgs_nmtf <- get_averages(avgs_nmtf, data_name, "results_nmtf", i, measures)
    avgs_nmtf_k_means <- get_averages(avgs_nmtf_k_means, data_name, "results_k_means_nmtf",i, measures)
}
#store results - one workbook for each test - each sheet a different measure
openxlsx::write.xlsx(avgs_nmtf, file = paste0(path_to_sim_folder, "/avgs_nmtf.xlsx"))
openxlsx::write.xlsx(avgs_nmtf_k_means, file = paste0(path_to_sim_folder, "/avgs_nmtf_k_means.xlsx"))
#list of files
#average col and row values, as well as total average
#now want to get averages and standard deviations across repetitions
#store results in a dataframe

results <- data.frame(Method=c(), Cluster=c(), Measure=c(), Mean=c(), Sd=c())

for (i in 1:length(measures)){
    results <- rbind(results,
                    avg_over_repeats(avgs_nmtf,i,measures[i],"SVD","Overall"),
                     avg_over_repeats(avgs_nmtf_k_means,i,measures[i],"k means","Overall"))
}
#change mean and sd columns to numeric 
colnames(results) <- c("Method", "Cluster", "Measure", "Mean", "S.d.")
results <-mutate_at(results,c("Mean","S.d."), as.numeric)
write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))

#make and saves plots

path_to_save <- paste0(paste0(path_to_sim_folder,"/"), paste0(paste0("Overall", "_diff_init"),".pdf"))
make_plot(results, "Overall", paste0("Overall", " performance with different initialisations"),
                paste0(n_repeats, " repetitions"),
                path_to_save)

path_to_save <- paste0(path_to_sim_folder, "/diff_init_table_")
#produce tables
all_table2(results, measures, path_to_save)
