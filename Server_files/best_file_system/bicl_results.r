args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
folder_list = args[2]
plot_title = args[3]
source("summary_funcs.r")
source(paste0(path_to_sim_folder, "/sim_parameters.r"))

#file to extract results
#read in folder_list
folder_names <- read.delim(folder_list,
         colClasses = "character", header = FALSE)[, 1]
n_repeats <- length(folder_names)
# n_repeats x 3 matrices
measures <- c("CSR", "Rel", "Rec", "F_score", "BiS")
averages <- c("Overall")
results <- data.frame(Method = c(), Choice = c(), Measure=c(), Mean=c(), Sd=c())
#analyse each method#
for (k in 1:length(file_names)){
    avgs <- replicate(length(measures),
                 matrix(nrow = n_repeats, ncol = length(averages)),
                  simplify = FALSE)
    for (i in 1:n_repeats){
        data_name <- paste0(paste0(path_to_sim_folder, "/data/"),
                folder_names[[i]])
        avgs <- get_averages(avgs,
                data_name, paste0(file_names[k], "_results"), i, measures)
    }
    for (i in 1:length(measures)){
        results <- rbind(results,
                         avg_over_repeats2(avgs, i,
                         measures[i], method_vec_res[k],
                          factor_vec[k], averages[1]))
    }
}
#change mean and sd columns to numeric
colnames(results) <- c("Method", factor_name,
                         "Choice", "Measure", "Mean", "S.d.")
results <- mutate_at(results, c("Mean", "S.d."), as.numeric)
results[["Choice"]] <- as.factor(results[["Choice"]])
results[["Method"]] <- as.factor(results[["Method"]])
results[[factor_name]] <- as.factor(results[[factor_name]])
measures <- c("CSR", "Rel", "Rec", "F_score", "BiS")
results[["Measure"]] <- factor(results[["Measure"]], levels = measures)

write.csv(results, paste0(path_to_sim_folder, "/all_results.csv"))

#now produce plots
#all measures, only ResNMTF
colnames(results) <- c("Method", "Factor", "Choice", "Measure", "Mean", "S.d.")
results_nmtf <- subset(results, Method == "ResNMTF")
measures <- c("CSR", "Relevance", "Recovery", "F score", "BiS")
if(phi_constant){
        subhead <- TeX(paste0(kept_factor,
                 "$\\phi$", " = ", val, ", ",
                  n_repeats, " repetitions."))
}else{
        subhead <- TeX(paste0(kept_factor,
                  n_repeats, " repetitions."))
}

path_to_save <- paste0(path_to_sim_folder, "/resNMTF_plot_", factor, ".pdf")
make_plot_line(results_nmtf, "Overall", plot_title,
                subhead, x_title, "Measure",
                path_to_save, measures)
if(phi_constant){
        #now produce one plot for each measure
        measures <- c("Relevance", "Recovery", "F score", "BiS", "CSR")
        n_measures <- length(measures)
        measure_vec <-  c("Rel", "Rec", "F_score", "BiS", "CSR")
        plots <- vector("list", length = n_measures)
        tables <- vector("list", length = n_measures)
        for(i in 1:n_measures){
                path_to_save_plots <-   paste0(path_to_sim_folder,
                         "/", measure_vec[i], "_plot.pdf")
                plots[[i]] <- make_plot_measure(results, "Overall", plot_title,
                                subhead, x_title, measures[i],
                                path_to_save_plots, measure_vec[i],
                                c("GFA", "iSSVD", "NMTF", "ResNMTF"))
                path_to_save_tables <- paste0(path_to_sim_folder,
                         "/", measure_vec[i], "_tbl")
                tables[[i]] <- measure_table(results, measure_vec[i], 
                        path_to_save_tables, n_repeats, col_names_tables,
                         subhead, x_title, ndig = 4)
        }
        #combine plots
        path_to_save_plot <- paste0(path_to_sim_folder, "/all_plot_", factor, ".pdf")
        all_plots(plots, path_to_save_plot, n_col_plots, n_row_plots)
        #combine tables
        path_to_save_table <- paste0(path_to_sim_folder, "/all_tables_", factor)
        combine_tables(tables, path_to_save_table, measures, subhead, x_title)
}

#produce correlation table
path_to_save <- paste0(path_to_sim_folder, "/", factor)
cor_res <- corr_table(results, path_to_save)

#produce table 1 - all measures for resNMTF
path_to_save <- paste0(path_to_sim_folder, "/resNMTF_table_", factor)
measure_labels <- c("F_score",  "Rel", "Rec", "CSR",  "BiS")
col_names <- c("Measure", "F score", "Relevance", "Recovery", "CSR", "BiS")
all_table2(results_nmtf, measure_labels, path_to_save,
                 n_repeats, col_names, subhead, table_head = x_title, feat = "Factor")
