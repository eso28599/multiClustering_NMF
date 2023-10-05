args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
folder_list = args[2]
plot_title = args[3]
source("summary_funcs.r")
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(kableExtra))
library(stringr)

#file to extract results
#read in folder_list
folder_names <- read.delim(folder_list,
         colClasses = "character", header = FALSE)[, 1]

n_repeats <- length(folder_names)
# n_repeats x 3 matrices
measures <- c("CSR", "Rel", "Rec", "F_score", "BiS", "Restrictions")
averages <- c("Overall")
results <- data.frame(Method = c(), Choice = c(), Measure=c(), Mean=c(), Sd=c())

method_vec <- paste0("res_nmtf_", phi_vals)
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
                         measures[i], 0.01*as.numeric(str_replace(method, "res_nmtf_", "")), averages[1]))
    }
}
#change mean and sd columns to numeric
colnames(results) <- c("Method", "Choice", "Measure", "Mean", "S.d.")
results <- mutate_at(results, c("Mean", "S.d."), as.numeric)
results[["Method"]] <- factor(results[["Method"]], levels = phi_vals*0.01)
results[["Choice"]] <- as.factor(results[["Choice"]])
measures <- c("CSR", "Rel", "Rec", "F_score", "BiS", "Restrictions")
#results[["Method"]] <- mutate_at(results, "Method"
                 #function(x) as.numeric(str_replace(x, "res_nmtf_", ""))*0.01)
#results[["Method"]] <- factor(results[["Method"]], levels = method_vec)
#for (i in 1:length(phi_vals)){
 #       results[["Method"]] <- ifelse(results[["Method"]] == method_vec[i],
  #              phi_vals[i], results[["Method"]])
        #results[["Method"]][results[["Method"]] == method_vec[i]] <- as.character(phi_vals[i]*0.01)
#}
#results[["Method"]] <- factor(results[["Method"]], levels = as.character(phi_vals * 0.01))
results2 <- subset(results, Measure %in% c("CSR", "F_score", "BiS", "Restrictions"))
results2[["Measure"]] <- factor(results2[["Measure"]], 
        levels = c("CSR", "F_score", "BiS", "Restrictions"))
write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))

#make and saves plot 

number <- length(row_cl_dims)
number_bicls <- length(row_start[[1]])
subhead <- paste0(paste0(number, " views, "), paste0(number_bicls, " biclusters, ") ,paste0(n_repeats, " repetitions."))
plot_title <- "The effect of increasing the strength of regularisation on performance"
x_title <- "Value of phi"
path_to_save <- paste0(path_to_sim_folder, "/increasing_phi_plot.pdf")
suppressMessages(make_plot_line_phi(results2, "Overall", plot_title,
                subhead,x_title,
                path_to_save))

path_to_save <- paste0(path_to_sim_folder, "/increasing_phi_")
#produce tables

col_names <- c("Phi", "CSR", "Relevance", "Recovery", "F score", "BiS", "Restrictions")
#subhead <- paste0(paste0(n_clusts, "biclusters, ") ,paste0(n_repeats, " repetitions."))
all_table2(results, measures, path_to_save, n_repeats, col_names, subhead)
