#dna microarray data
source("resNMTF_funcs.r")
micro_dna <- read.csv("pnas2.csv", header = T)
phi_mat <- matrix(0, 1, 1)
apply_data <-  restMultiNMTF_run(Xinput = list(t(micro_dna)), k_min=3,
            phi = phi_mat, relErr = "Sil")
colSums(apply_data$row_clusters[[1]])
colSums(apply_data$col_clusters[[1]])
(1:156)[apply_data$row_clusters[[1]][,2]==1]
length(intersect(
    (intersect((1:156)[apply_data$col_clusters[[1]][,2]==1],
    (1:156)[apply_data$col_clusters[[1]][,3]==1])),(1:156)[apply_data$col_clusters[[1]][,1]==1]))
length(intersect((1:675)[apply_data$col_clusters[[1]][,2]==1],
    (1:675)[apply_data$col_clusters[[1]][,3]==1]))
apply(apply_data$Foutput[[1]]>(1/156), 2, sum)
apply_data2 <-  restMultiNMTF_run(Xinput = list(t(micro_dna)), k_min=2,
            phi = phi_mat, relErr = "Sil")
colSums(apply_data2$row_clusters[[1]])
colSums(apply_data2$col_clusters[[1]])
apply_data3 <-  restMultiNMTF_run(Xinput = list(t(micro_dna)), KK = c(2), LL = c(2),
            phi = phi_mat, relErr = "Sil")
(1:156)[apply_data3$row_clusters[[1]][,2]==1]
(1:156)[apply_data3$Foutput[[1]][,1]>1/156]

#consider only the AD samples
micro_dna2 <- t(read.csv("pnas2.csv", header = T))[1:(156-17), ]
apply_data <-  restMultiNMTF_run(Xinput = list(micro_dna2), k_min=3,
            phi = phi_mat, relErr = "Sil")
colSums(apply_data$row_clusters[[1]])
colSums(apply_data$col_clusters[[1]])
apply(apply_data$Foutput[[1]]>(1/139), 2, sum)
(1:139)[apply_data$Foutput[[1]][,5]>1/139]
micro_dna2 <- single_alt_l1_normalisation(micro_dna2)$newMatrix
sd_vec <- apply(micro_dna2, 2 , sd)
micro_dna2 <- micro_dna2[, sd_vec>0.01]
apply_data <-  restMultiNMTF_run(Xinput = list(micro_dna2), k_min=2,
            phi = phi_mat, relErr = "Sil")
colSums(apply_data$Foutput[[1]]>1/139)
(1:139)[(apply_data$Foutput[[1]][,1]>1/139)]


python_data <- import_matrix("results/data_issvd.xlsx")
python_rows <- import_matrix("results/row_issvd.xlsx")
python_cols <- import_matrix("results/col_issvd.xlsx")
phi_mat2 <- matrix(0, 2, 2)
phi_mat2[1,2] <- 200
new_data <- list("v1" = abs(python_data[[1]]), "v2" =abs(python_data[[2]]))
python_res <- restMultiNMTF_run(Xinput = python_data, k_min=3,
            phi = phi_mat2, relErr = "Sil")
bicl1_v1 <- (1:200)[python_res$row_clusters[[1]][,1]==1]
bicl2_v1 <- (1:200)[python_res$row_clusters[[1]][,2]==1]
bicl3_v1 <- (1:200)[python_res$row_clusters[[1]][,3]==1]
bicl4_v1 <- (1:200)[python_res$row_clusters[[1]][,4]==1]

bicl1_v2 <- (1:200)[python_res$row_clusters[[2]][,1]==1]
bicl2_v2 <- (1:200)[python_res$row_clusters[[2]][,2]==1]
bicl3_v2 <- (1:200)[python_res$row_clusters[[2]][,3]==1]
bicl4_v2 <- (1:200)[python_res$row_clusters[[2]][,4]==1]
source("evaluation_funcs.r")
python_res$row_clusters
evaluate_simulation_comp(python_res$row_clusters, python_res$col_clusters,
        python_rows, python_cols, python_data)

#bicluster one, rows 
r1 <- (1:200)[python_rows[[1]][,1]==1]
intersect(r1, bicl4_v1)

#attempt with positive data
python_res2 <- restMultiNMTF_run(Xinput = new_data, k_min=3,
            phi = phi_mat2, relErr = "Sil")
colSums(python_res2$row_clusters[[1]])
evaluate_simulation_comp(python_res2$row_clusters, python_res2$col_clusters,
        python_rows, python_cols, python_data)


python_res100 <- restMultiNMTF_run(Xinput = new_data, k_min=3,
            phi = 0.5*phi_mat2, relErr = "Sil")
colSums(python_res2$row_clusters[[1]])
evaluate_simulation_comp(python_res100$row_clusters, python_res100$col_clusters,
        python_rows, python_cols, python_data)


python_res0 <- restMultiNMTF_run(Xinput = new_data, k_min=3,
            phi = 0*phi_mat2, relErr = "Sil")
colSums(python_res2$row_clusters[[1]])
evaluate_simulation_comp(python_res0$row_clusters, python_res0$col_clusters,
        python_rows, python_cols, python_data)


python_res300 <- restMultiNMTF_run(Xinput = new_data, k_min=3,
            phi = 2/3*phi_mat2, relErr = "Sil")
colSums(python_res2$row_clusters[[1]])
evaluate_simulation_comp(python_res300$row_clusters, python_res300$col_clusters,
        python_rows, python_cols, python_data)

phi_val <- c(0, 10, 20, 30, 50, 75, 100, 200, 300, 500, 1000, 5000, 10000)
phi_mat2 <- matrix(0, 2, 2)
phi_mat2[1,2] <- 1
python_res <- restMultiNMTF_run(Xinput = new_data, k_min=3,
            phi = phi_val[1]*phi_mat2, relErr = "Sil") 
results <- vector("list", length=length(phi_val))
for(i in 1:length(phi_val)){
    python_res <- restMultiNMTF_run(Xinput = new_data, k_min=3,
            phi = phi_val[i]*phi_mat2, relErr = "Sil") 
    results[[i]] <-  evaluate_simulation_comp(python_res$row_clusters, python_res$col_clusters,
        python_rows, python_cols, python_data)   
}
f_score <- rep(0, 9)
biS <- rep(0,9)

for(i in 1:9){
    f_score[i] <- mean(results[[i]]$F_score)
    biS[i] <- mean(results[[i]]$BiS)
}

#load results
issvd_rows <- import_matrix("results/row_res_issvd.xlsx")
issvd_cols <- import_matrix("results/col_res_issvd.xlsx")
issvd_scores <- evaluate_simulation_comp(issvd_rows, issvd_cols,
        python_rows, python_cols, python_data)

issvd_rows_og <- import_matrix("results/row_res_issvd.xlsx")
issvd_cols_og <- import_matrix("results/col_res_issvd.xlsx")
issvd_scores_og <- evaluate_simulation_comp(issvd_rows_og, issvd_cols_og,
        python_rows, python_cols, python_data)

source("gfa_funcs.r")
gfa_res <- gfa_apply(new_data, 9)
gfa_res2 <- gfa_apply(new_data, 4)
gfa_scores2 <- evaluate_simulation_comp(gfa_res2$row_clusters, gfa_res2$col_clusters,
        python_rows, python_cols, python_data)
gfa_scores <- evaluate_simulation_comp(gfa_res$row_clusters, gfa_res$col_clusters,
        python_rows, python_cols, python_data)
#try GFA with no informative prior 
methods <- c("ResNMTF", "NMTF", "GFA", "iSSVD", "iSSVD_abs")
f_s <- c(f_score[7], f_score[1], mean(gfa_scores$F_score), 
        mean(issvd_scores_og$F_score), mean(issvd_scores$F_score) )
b_vec <- c(biS[7], biS[1], mean(gfa_scores$BiS),  mean(issvd_scores_og$BiS), mean(issvd_scores$BiS))
other_data_results <- data.frame("Method"=methods, 
        "F_score"=f_s, "BiS" = b_vec)
