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
#measures <- c("Accuracy", "NMI", "ARI")
measures <- c("CSR", "Rel", "Rev", "F_score", "BiS")
averages <- c("Overall")
results <- data.frame(Method = c(), Choice = c(), Measure=c(), Mean=c(), Sd=c())

#method_vec <- c("res_nmtf", "mix_nmtf", "issvd")
method_vec <- c("res_nmtf", "issvd")
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
#results[["Measure"]] <- ifelse(results[["Measure"]] == "Accuracy",
                    #"CSR", results[["Measure"]])
#measures <- c("CSR", "NMI", "ARI")
measures <- c("CSR", "Rel", "Rev", "F_score", "BiS")
results[["Measure"]] <- factor(results[["Measure"]], levels = measures)
results[["Method"]] <- ifelse(results[["Method"]] == "res_nmtf",
                    "ResNMTF", ifelse(results[["Method"]] == "mix_nmtf",
                    "ResNMTF (RelErr)", "iSSVD"))
write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))

#make and saves plots
#change to phi plot, giving lines for different methods
path_to_save <- paste0(path_to_sim_folder, "/scen.pdf")
suppressMessages(make_plot(results, "Overall", plot_title,
                paste0(n_repeats, " repetitions"),
                path_to_save))

path_to_save <- paste0(path_to_sim_folder, "/scenario_")
#produce tables
all_table2(results, measures, path_to_save, n_repeats, group = FALSE)
