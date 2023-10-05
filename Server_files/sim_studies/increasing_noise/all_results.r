args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
folder_list = args[2]
plot_title = args[3]
source("summary_funcs.r")
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(kableExtra))

#file to extract results
#read in folder_list
folder_names <- read.delim(folder_list,
         colClasses = "character", header = FALSE)[, 1]

n_repeats <- length(folder_names)
# n_repeats x 3 matrices
measures <- c("CSR", "Rel", "Rec", "F_score", "BiS")
averages <- c("Overall")
results <- data.frame(Method = c(), Choice = c(), Measure=c(), Mean=c(), Sd=c())

method_vec <- paste0("res_nmtf_", 1:10)
#method_vec <- c("res_nmtf_2", "res_nmtf_3", "res_nmtf_4", "res_nmtf_5")
#analyse each method#
for (method in method_vec){
    avgs <- replicate(length(measures),
                 matrix(nrow = n_repeats, ncol = length(averages)),
                  simplify = FALSE)
    for (i in 1:n_repeats){
        data_name <- paste0(paste0(path_to_sim_folder, "/data/"),
                folder_names[[i]])
        avgs <- get_averages(avgs,
                data_name, paste0(method, "_results"), i, measures)
    }
    for (i in 1:length(measures)){
        results <- rbind(results,
                         avg_over_repeats(avgs, i,
                         measures[i], method, averages[1]))
    }
}
#change mean and sd columns to numeric
colnames(results) <- c("Method", "Choice", "Measure", "Mean", "S.d.")
results <- mutate_at(results, c("Mean", "S.d."), as.numeric)
results[["Method"]] <- factor(results[["Method"]], levels = method_vec)
results[["Choice"]] <- as.factor(results[["Choice"]])
measures <- c("CSR", "Rel", "Rec", "F_score", "BiS")
for (i in 1:10){
        results[["Method"]] <- ifelse(results[["Method"]] == method_vec[i],
         i, results[["Method"]])
}
results2 <- subset(results, Measure %in% c("CSR", "F_score", "BiS"))
results2[["Measure"]] <- factor(results2[["Measure"]], 
        levels = c("CSR", "F_score", "BiS"))
write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))

#make and saves plot 
#
plot_title <- "The effect of increasing the level of noise on performance"
path_to_save <- paste0(path_to_sim_folder, "/increasing_noise_plot.pdf")
suppressMessages(make_plot_line(results2, "Overall", plot_title,
                paste0(n_repeats, " repetitions"), "Noise level",
                path_to_save))

path_to_save <- paste0(path_to_sim_folder, "/increasing_noise")
#produce tables
all_table2(results, measures, path_to_save, n_repeats, group = FALSE)
