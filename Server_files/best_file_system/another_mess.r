#increasing views
library("readxl")
library("rio")
source("data_generation.r")
source("resNMTF_funcs.r")
source("evaluation_funcs.r")
import_matrix <- function(filename, col_n = TRUE){
    return(lapply(import_list(filename, col_names=col_n), function(x) as.matrix(x)))
}
# 3 views, 3 biclusters
source("best_file_system/increasing_views/sim_parameters.r")
data_all <- multi_view(row_cl_dims, col_cl_dims, row_start, row_end,  col_start, col_end,
        path_to_data_folder, row_same_shuffle, col_same_shuffle, noise_level)
data_views <- data_all$data_views

#10
t1 <- Sys.time()
nmtf_results <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]]/50000, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results$row_clusters[[i]],nmtf_results$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}
#50000
t1 <- Sys.time()
nmtf_results_500000 <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]], relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_500000$row_clusters[[i]],nmtf_results_500000$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}

#1
t1 <- Sys.time()
nmtf_results_1 <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]]/500000, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_1$row_clusters[[i]],nmtf_results_1$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}

#now using original version of star product

#10
t1 <- Sys.time()
nmtf_results_o <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]]/50000, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_o$row_clusters[[i]],nmtf_results_o$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}
#50000
t1 <- Sys.time()
nmtf_results_500000_o <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]], relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_500000_o$row_clusters[[i]],nmtf_results_500000_o$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}

#1
t1 <- Sys.time()
nmtf_results_1_o <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]]/500000, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_1_o$row_clusters[[i]],nmtf_results_1_o$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}

#with subsetted sil score
#50000
t1 <- Sys.time()
nmtf_results_500000_s <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]], relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_500000_s$row_clusters[[i]],nmtf_results_500000_s$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}

#50000
t1 <- Sys.time()
nmtf_results_500000000_s <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]]*1000, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_500000000_s$row_clusters[[i]],nmtf_results_500000000_s$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
    print(colSums(nmtf_results_500000000_s$row_clusters[[i]]))
    print(colSums(nmtf_results_500000000_s$col_clusters[[i]]))
}


#all udpates
t1 <- Sys.time()
nmtf_results_up <- restMultiNMTF_run(Xinput = data_views[1:3],
         phi = phi_mat[2][[1]]/500000, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
for(i in 1:3){
    print(jaccard_results(nmtf_results_up$row_clusters[[i]],nmtf_results_up$col_clusters[[i]],
         data_all$truth_rows[[i]], data_all$truth_cols[[i]]))
}
evaluate_simulation_comp(nmtf_results_up$row_clusters, nmtf_results_up$col_clusters,
 data_all$truth_rows[1:3],  data_all$truth_cols[1:3], data_views[1:3], index = 2)
sil_score(data_views[[3]], nmtf_results_up$row_clusters[[3]],
                                  nmtf_results_up$col_clusters[[3]], 2)


#look at the non-working for more than 6 
data_views <- import_matrix("results_24_08_p200/increasing_bicl2/data.xlsx", col_n = TRUE)
rows <- import_matrix("results_24_08_p200/increasing_bicl2/true_rows.xlsx", col_n = TRUE)
cols <- import_matrix("results_24_08_p200/increasing_bicl2/true_cols.xlsx", col_n = TRUE)
image(data_views[[1]])
reorder_rows <- 
image(data_views[[1]])
#nmtf
phi <- matrix(0, nrow = 3, ncol = 3)
phi[1, c(2,3)] <- 1
phi[2, c(3)] <- 1
phi_mat <- 200 * phi
nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat, relErr = "Sil")
evaluate_simulation_comp(nmtf_results$row_clusters, nmtf_results$col_clusters,
 rows,  cols, data_views, index = 2)


data_views2 <- import_matrix("results_24_08_p200/increasing_bicl4/data.xlsx", col_n = TRUE)
rows2 <- import_matrix("results_24_08_p200/increasing_bicl4/true_rows.xlsx", col_n = TRUE)
cols2 <- import_matrix("results_24_08_p200/increasing_bicl4/true_cols.xlsx", col_n = TRUE)
image(data_views[[1]])
reorder_rows <- 
image(data_views[[1]])
#nmtf
phi <- matrix(0, nrow = 3, ncol = 3)
phi[1, c(2,3)] <- 1
phi[2, c(3)] <- 1
phi_mat <- 200 * phi
nmtf_results2 <- restMultiNMTF_run(Xinput = data_views2, phi = phi_mat, relErr = "Sil")
evaluate_simulation_comp(nmtf_results2$row_clusters, nmtf_results2$col_clusters,
 rows2,  cols2, data_views2, index = 2)


row_cl_dims <- rep(200, 4)

#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250, 150)
noise_level <-5
row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 3 biclusters
row_start <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 200), c(75, 150, 200),
            c(75, 150, 200), c(75, 150, 200))
col_start <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151), c(1, 51, 101))
col_end <- list(c(30, 60, 100), c(10, 30, 50),
                 c(100, 150, 250), c(50, 100, 150))
data_all <- multi_view(row_cl_dims, col_cl_dims, row_start, row_end,  col_start, col_end,
        path_to_data_folder, row_same_shuffle, col_same_shuffle, noise_level)
data_views <- data_all$data_views
#nmtf
phi <- matrix(0, nrow = 4, ncol = 4)
phi[1, c(2, 3, 4)] <- 1
phi[2, c(3, 4)] <- 1
phi[3, 4] <- 1
phi_mat <- 20000 * phi
nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat, relErr = "Sil")
source("best_file_system/evaluation_funcs.r")
nmtf_scores2 <- evaluate_simulation_comp(nmtf_results$row_clusters , nmtf_results$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)

clust_results <- clustering_res_NMTF(data_views, nmtf_results$Foutput,
               nmtf_results$Goutput, nmtf_results$Soutput,
                3, 2, phi_mat, epsilon = 0.02)
colSums(apply(nmtf_results$Goutput[[i]],
          2, function(x) as.numeric(x > (1 / dim(nmtf_results$Goutput[[i]])[1]))))
nmtf_results_0 <- restMultiNMTF_run(Xinput = data_views, phi = 0*phi_mat, relErr = "Sil")

#2 biclusters
row_start <- list(c(1, 76), c(1, 76), c(1, 76),c(1, 76))
row_end <- list(c(75, 200), c(75, 200),
            c(75, 200),c(75, 200))
col_start <- list(c(1, 31), c(1, 11), c(1, 101), c(1, 51))
col_end <- list(c(30, 100), c(10, 50),
                 c(100, 250), c(50,  150))
data_all2 <- multi_view(row_cl_dims, col_cl_dims, row_start, row_end,  col_start, col_end,
        path_to_data_folder, row_same_shuffle, col_same_shuffle, noise_level)
data_views2 <- data_all2$data_views
nmtf_results2 <- restMultiNMTF_run(Xinput = data_views2, phi = phi_mat, relErr = "Sil")

example_data <- one_view_adv(200, 100, c(1, 76, 151), c(75, 150, 200),  c(1, 31, 61), c(30, 60, 100), overlap = FALSE, noise=5)
example_data_shuffled <- example_data$view[sample(200), sample(100)]
image(t(example_data_shuffled))
image(t(example_data$view))

bicl_res <- read.csv("best_file_system/increasing_bicl/all_results.csv")
#make and saves plot
source("best_file_system/increasing_bicl/sim_parameters.r")
n_repeats <- 15
number <- length(row_cl_dims)
subhead <- TeX(paste0(number, " views, ",
                 "$\\phi$", " = ", val * 0.01, ", ",
                  n_repeats, " repetitions."))
x_title <- "Number of biclusters"
plot_title <- "The effect of increasing the number of biclusters on performance"
path_to_save <- "test.pdf"
suppressMessages(make_plot_line(bicl_res, "Overall", plot_title,
                subhead, x_title,
                path_to_save))

path_to_save <- paste0(path_to_sim_folder, "/increasing_biclusts")

file_path <- "best_file_system/increasing_phi4/res_nmtf_42000"
row_filename <- paste0(file_path, "_row_clusts.xlsx")
col_filename <- paste0(file_path, "_col_clusts.xlsx")
f_sd_filename <- paste0(file_path, "_f_diff.csv")
    #results for nmtf
data_name <- "best_file_system/increasing_phi4"
true_rows <- import_matrix(paste0(data_name, "/true_rows.xlsx"))
true_cols <- import_matrix(paste0(data_name, "/true_cols.xlsx"))
data_views <- import_matrix(paste0(data_name, "/data.xlsx"), col_n = TRUE)
results <- evaluate_simulation_comp(import_matrix(row_filename),
                                              import_matrix(col_filename),
                                              true_rows,
                                               true_cols,
                                                data_views)
tol <-1e-06
f_sd_mat <- matrix(as.numeric(read.csv(f_sd_filename)$x < tol),
                 2, dim(results$CSR)[2])
    rownames(f_sd_mat) <- c("Row-clustering", "Column-clustering")
    results[["Restrictions"]] <- f_sd_mat

#c


normX <- normalizeData(data_all$data_views, type = "center")
opts <- getDefaultOpts(bicluster = TRUE)
opts <- getDefaultOpts()
gfa_res <- gfa(data_all$data_views,K =10, opts)
gfa_res <- gfa(normX[[1]],K =10, opts)
gfa_res$groups
gfa_res$X
gfa_res$Z
colSums(gfa_res$X>0)
colSums(gfa_res$Z[1:100,]>0)
test <-rep(list(apply(gfa_res$X>0, 2, as.numeric)), 4)
results <- gfa_apply(data_all$data_views)

source("evaluation_funcs.r")
res_gfa <- evaluate_simulation_comp(results$row_clusts, results$col_clusts,    
    data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)