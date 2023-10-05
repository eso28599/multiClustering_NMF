new_sil_score <- function(Xinput, Foutput, Goutput, row_clustering, col_clustering, index){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    n_views <- length(Xinput)
    n_clusts <- dim(Foutput[[1]])[2]
    relations <- matrix(0, nrow = n_views, n_clusts)
    sil_score <- matrix(0, nrow = n_views, n_clusts)
    if (index == 2) {
        clust_one <- col_clustering
        clust_two <- row_clustering
    }else{
        clust_one <- row_clustering
        clust_two <- col_clustering
    }
    for (i in 1:n_views){
        s_mat <- matrix(0, nrow = n_clusts, n_clusts)
        a_mat <- matrix(0, nrow = n_clusts, n_clusts)
        b_mat <- matrix(0, nrow = n_clusts, n_clusts)
        for (k in 1:n_clusts){
          #select data from specific column clustering
          if (sum(clust_one[[i]][, k] == 1) == 0){
            #print("no cols")
            s_mat[k, ] <- 0
          }else{
            if (index == 2){
              new_data <- Xinput[[i]][, (clust_one[[i]][, k] == 1)]
            }else{
              new_data <- t(Xinput[[i]][(clust_one[[i]][, k] == 1), ])
            }
            spear_dists <- as.matrix(dist(new_data, upper = FALSE, diag = TRUE))/ (dim(new_data)[2])
            for (j in 1:n_clusts){
              indices <- clust_two[[i]][, j] == 1
              if (sum(indices) == 0) {
                #("none")
                s_mat[k, j] <- 0
              }else{
                #a_vals <- rowMeans(spear_dists[indices, indices])
                a_vals <- apply(spear_dists[indices, indices], 1, function(x) sum(x)/(length(x)-1))
                #other clusts
                other <- (1:n_clusts)[-j]
                b_vec <- c()
                b_vals <- vector("list", length = (n_clusts - 1))
                t <- 1
                for(l in other){
                    oth_ind <- clust_two[[i]][, l] == 1
                    b_val <- rowMeans(spear_dists[indices, oth_ind])
                    b_vec <- c(b_vec, mean(b_val))
                    b_vals[[t]] <- b_val
                    t <- t + 1
                }
                closest <- which.min(b_vec)
                print(b_vec)
                b_vals <- b_vals[[closest]]
                s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
                s_mat[k, j] <- mean(s_con)
                #a_mat[k, j] <- mean(a_vals)
                #b_mat[k,j] <- mean(b_vals)
               }
            }
          }
        }
        print(s_mat)
        if(n_clusts == 2){
          relations[i, ] <- apply(s_mat, 1, which.min)
          sil_score[i, ] <- apply(s_mat, 1, min)
        }else{
          relations[i, ] <- apply(s_mat, 1, which.max)
          sil_score[i, ] <- apply(s_mat, 1, max)
        }
        
      
      #find the contribution to the s_score of the biclusters from this view
      
    }
  acc_vals <- sil_score[sil_score!=0]
  sil <- ifelse(n_clusts == 1, mean(acc_vals),
                 mean(acc_vals) - 2 * sd(acc_vals))
  #sil <- ifelse(n_clusts == 1, mean(sil_score[sil_score!=0]),
                 #mean(sil_score) - 2 * sd(sil_score))
  #return relationships and mean of sil_score across views
  return(list("sil" = sil, "scores" = sil_score, "relations" = relations))
}


# 5 biclusters
row_start_5 <- list(c(1, 51, 76, 111, 151), c(1, 51, 76, 111, 151),
                     c(1, 51, 76, 111, 151), c(1, 51, 76, 111, 151))
row_end_5 <- list(c(50, 75, 110, 150, 200), c(50, 75, 110, 150, 200),
                 c(50, 75, 110, 150, 200), c(50, 75, 110, 150, 200))
col_start_5 <- list(c(1, 21, 31, 46, 61), c(1, 6, 11, 21, 31),
                         c(1, 51, 101, 126, 151), c(1, 26, 51, 76, 101))
col_end_5 <- list(c(20, 30, 45, 60, 100), c(5, 10, 20, 30, 50),
                     c(50, 100, 125, 150, 250), c(25, 50, 75, 100, 150))
noise_level <- 5
#each dataset with 200 samples, same clusters across views
row_cl_dims <- rep(200, 4)

#100, 50,250, 150 features respectively
col_cl_dims <- c(100, 50, 250, 150)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE
niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_cl_dims)
no_row_cl <- sapply(row_cl_dims, length)
no_col_cl <- sapply(col_cl_dims, length)
#parameters for restrictive NMTF
val <- 200
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, c(2, 3, 4)] <- 1
phi[2, c(3, 4)] <- 1
phi[3, 4] <- 1
phi_mat <- val * phi

data_all <- multi_view(row_cl_dims,
        col_cl_dims, row_start_5, row_end_5, col_start_5, col_end_5, 
        noise=noise_level, row_same_shuffle=TRUE, col_same_shuffle = FALSE, seed=10)
source("best_file_system/resNMTF_funcs.r")
nmtf_results <- restMultiNMTF_run(Xinput = data_all$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")
#evaluate results
source("evaluation_funcs.r")
res0 <- evaluate_simulation_comp(nmtf_results$row_clusters , nmtf_results$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)
source("resNMTF_funcs.r")
orignal_method <- restMultiNMTF_run(Xinput = data_all$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")
source("evaluation_funcs.r")
res1 <- evaluate_simulation_comp(orignal_method$row_clusters , orignal_method$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)

# 4 biclusters
row_start_4 <- list(c(1, 76, 111, 151), c(1, 76, 111, 151),
                     c(1, 76, 111, 151), c(1, 76, 111, 151))
row_end_4 <- list(c(75, 110, 150, 200), c(75, 110, 150, 200),
                 c(75, 110, 150, 200), c(75, 110, 150, 200))
col_start_4 <- list(c(1, 31, 46, 61), c(1, 11, 21, 31),
                 c(1, 101, 126, 151), c(1, 51, 76, 101))
col_end_4 <- list(c(30, 45, 60, 100), c(10, 20, 30, 50),
             c(100, 125, 150, 250), c(50, 75, 100, 150))

jaccard_results(new_method_4$row_clusters[[3]][,1:2], new_method_4$col_clusters[[3]][, 1:2],
 data_all$truth_rows[[3]], data_all$truth_cols[[3]])



rows <- list()
rows[[1]]
source("resNMTF_funcs.r")
data_all <- multi_view(row_cl_dims,
        col_cl_dims, row_start_4, row_end_4, col_start_4, col_end_4, 
        noise=noise_level, row_same_shuffle=TRUE, col_same_shuffle = FALSE, seed=10)
nmtf_results4 <- restMultiNMTF_run(Xinput = data_all$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")

nmtf_results4 <- restMultiNMTF_run(Xinput = list(data_all$data_views[[1]]), k_min=3,  phi=  matrix(0,1,1), relErr = "Sil")

single_nmtf <- single_nmtf(data_all$data_views)
source("evaluation_funcs.r")
res_single <- evaluate_simulation_comp(single_nmtf$row_clusters , single_nmtf$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)
#evaluate results
source("evaluation_funcs.r")
res_org_4 <- evaluate_simulation_comp(nmtf_results4$row_clusters , nmtf_results4$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)
source("resNMTF_funcs.r")
new_method_4 <- restMultiNMTF_run(Xinput = data_all$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")
source("evaluation_funcs.r")
res1 <- evaluate_simulation_comp(new_method_4$row_clusters , new_method_4$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)
#test with lower stab threshold
source("resNMTF_funcs.r")
new_method_4_lt <- restMultiNMTF_run(Xinput = data_all$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")

#check with a lower phi
source("resNMTF_funcs.r")
new_method_4low <- restMultiNMTF_run(Xinput = data_all$data_views, k_min=3,
            phi = 10/200*phi_mat, relErr = "Sil")
source("evaluation_funcs.r")
res1_low <- evaluate_simulation_comp(new_method_4low$row_clusters , new_method_4low$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)
res1_low



X1 <- single_alt_l1_normalisation(data_all$data_views[[1]])$newMatrix
X1_hat <- nmtf_results4$Foutput[[1]]%*%nmtf_results4$Soutput[[1]]%*%t(nmtf_results4$Goutput[[1]])

row_start_2 <- list(c(1, 76), c(1, 76), c(1, 76),  c(1, 76))
row_end_2 <- list(c(75, 200), c(75, 200), c(75, 150), c(75, 150))
col_start_2 <- list(c(1, 31), c(1, 11), c(1, 101), c(1, 51))
col_end_2 <- list(c(30, 100), c(10, 50), c(100, 250), c(50, 150))


data2 <- multi_view(row_cl_dims,
        col_cl_dims, row_start_2, row_end_2, col_start_2, col_end_2, 
        noise=noise_level, row_same_shuffle=TRUE, col_same_shuffle = FALSE, seed=10)
source("resNMTF_funcs.r")
nmtf_results2 <- restMultiNMTF_run(Xinput = data2$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")


# 3 biclusters
row_start_3 <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end_3 <- list(c(75, 150, 200), c(75, 150, 200),
            c(75, 150, 200), c(75, 150, 200))
col_start_3 <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151), c(1, 51, 101))
col_end_3 <- list(c(30, 60, 100), c(10, 30, 50),
                 c(100, 150, 250), c(50, 100, 150))
data3 <- multi_view(row_cl_dims,
        col_cl_dims, row_start_3, row_end_3, col_start_3, col_end_3, 
        noise=noise_level, row_same_shuffle=TRUE, col_same_shuffle = FALSE, seed=10)
nmtf_results3 <- restMultiNMTF_run(Xinput = data3$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")


nmtf_results_4_low <- restMultiNMTF_run(Xinput = data_all$data_views, k_min=3,
            phi = phi_mat, relErr = "Sil")
