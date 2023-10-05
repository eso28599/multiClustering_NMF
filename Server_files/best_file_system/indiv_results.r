args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
folder_list = args[2]
plot_title = args[3]
source("summary_funcs.r")
source(paste0(path_to_sim_folder, "/sim_parameters.r"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(latex2exp))

#file to extract results
#read in folder_list
folder_names <- read.delim(folder_list,
         colClasses = "character", header = FALSE)[, 1]

n_repeats <- length(folder_names)
# n_repeats x 3 matrices
measures <- c("CSR", "Rel", "Rec", "F_score", "BiS")
averages <- c("Overall")
results <- data.frame(Method = c(), Choice = c(), Measure=c(), Mean=c(), Sd=c())

method_vec <- c("res_nmtf_50", "res_nmtf_200", "res_nmtf_300", "res_nmtf_500", "res_nmtf_1000")
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
results[["Method"]] <- ifelse(results[["Method"]] == "res_nmtf_50",
                    "50", ifelse(results[["Method"]] == "res_nmtf_200",
                    "200",  ifelse(results[["Method"]] == "res_nmtf_300",
                    "300",  ifelse(results[["Method"]] == "res_nmtf_500",
                    "500",  "1000"))))
results[["Method"]] <- factor(results[["Method"]],
                         levels = c("50", "200", "300", "500", "1000"))              
results2 <- subset(results, Measure %in% c("CSR", "F_score", "BiS"))
results2[["Measure"]] <- factor(results2[["Measure"]], 
        levels = c("CSR", "F_score", "BiS"))
write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))

#make and saves plot 
number <- length(row_cl_dims[[1]])
number_bicls <- length(row_start_50[[1]])
subhead <- TeX(paste0(number, " views, ", number_bicls, " biclusters, ",
                 "$\\phi$", " = ", val, ", ",
                  n_repeats, " repetitions."))
x_title <- "Number of individuals"
plot_title <- "The effect of increasing the number of individuals on performance"
path_to_save <- paste0(path_to_sim_folder, "/increasing_indiv_plot.pdf")
suppressMessages(make_plot_line(results2, "Overall", plot_title,
         subhead, x_title,
                path_to_save))

path_to_save <- paste0(path_to_sim_folder, "/increasing_indiv")
col_names <- c("No. of individuals","CSR", "Relevance","Recovery", "F score", "BiS")
#produce tables
all_table2(results, measures, path_to_save, n_repeats, col_names, subhead, group = FALSE)
