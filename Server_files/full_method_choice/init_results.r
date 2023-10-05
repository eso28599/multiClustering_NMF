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
end_vec <- c("Iter", "Conv")
choice_vec <- c("Error", "RelErr", "Sil")

measures <- c("accuracy", "NMI", "ARI")
averages <- c("Overall")
results <- data.frame(Method=c(), Choice=c(), Measure=c(), Mean=c(), Sd=c())

for (end in end_vec){
    avgs_nmtf <- replicate(length(choice_vec), matrix(nrow = n_repeats, ncol = length(averages)), simplify=FALSE)
    for (choice in choice_vec){
        for (i in 1:length(folder_names)){
            data_name <- paste0(paste0(path_to_sim_folder, "/data/"), folder_names[[i]])
            avgs_nmtf <- get_averages(avgs_nmtf, data_name, paste0(paste0(end, choice),"_results_nmtf"), i, measures)
        }
        for (i in 1:length(measures)){
            results <- rbind(results,
                            avg_over_repeats(avgs_nmtf,i,measures[i],end,choice))
            }
    }
}

#change mean and sd columns to numeric 
colnames(results) <- c("Method", "Choice", "Measure", "Mean", "S.d.")
results <-mutate_at(results,c("Mean","S.d."), as.numeric)
results <-mutate_at(results,c("Method", "Choice"), as.factor)
write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))


#make and saves plots
##checked up to here!!! need to add make_plot_params to summary_funcs file and change so works 
#without subsetting by cluster
for(measure in measures){
    path_to_save <- paste0(paste0(path_to_sim_folder,"/"), paste0(paste0(measure, "_diff_selec"),".pdf"))
    make_plot_params(results, measure, paste0(measure, " with different selection methods"),
                paste0(n_repeats, " repetitions"),
                path_to_save)
}

path_to_save <- paste0(path_to_sim_folder, "/diff_init_table_")
#produce tables
all_table2(results, measures, path_to_save, group=TRUE)
