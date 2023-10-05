suppressPackageStartupMessages(library(mclust))
library(clue)
library(aricode)
library(cluster)
suppressPackageStartupMessages(library(proxy))

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <-  length(a) + length(b) - intersection
    if (union == 0 ){
      return(0)
    }else {
      return(intersection / union)
    }
}
cart_prod <- function(a, b) {
  prod <- c()
  # check a or b are not empty sets
  if (length(a) == 0 || length(b) == 0) {
    return(NULL)
  }else{
    for (k in 1:length(a)){
      prod <- c(prod, paste(a[k], b))
    }
  return(prod)
  }
}
jaccard_results <- function(row_c, col_c, true_r, true_c, stability = FALSE){
  m <- sum(colSums(row_c) != 0) #no of clusters detected
  n <- sum(colSums(true_r) != 0) #no of true clusters
  samps <- 1:nrow(row_c)
  feats <- 1:nrow(col_c)
  if ((m == 0 && n != 0) || (n == 0 && m != 0)) {
    return(list("rec" = rep(0, 2), "rel" = rep(0, 2), "f_score" = rep(0, 2)))
  }
  # if no biclusters present and none detected - score of 1
  if (m == 0 && n == 0) {
    return(list("rec" = rep(1, 2), "rel" = rep(1, 2), "f_score" = rep(1, 2)))
  }
  #initialise storage of jaccard index between pairs
  jac_mat <- matrix(0, nrow = m, ncol = n)
  for (i in 1:m){
    r_i <- samps[row_c[, i] == 1]
    c_i <- feats[col_c[, i] == 1]
    m_i <- cart_prod(r_i, c_i)
    for (j in 1:n){
        tr_i <- samps[true_r[, j] == 1]
        tc_i <- feats[true_c[, j] == 1]
        m_j <- cart_prod(tr_i, tc_i)
        jac_mat[i, j] <- jaccard(m_i, m_j)
    }
  }
  rel <- ifelse(sum(apply(jac_mat, 1, max)!=0)==0, 0, 
    sum(apply(jac_mat, 1, max))/sum(apply(jac_mat, 1, max)!=0))
  rec <- ifelse(sum(apply(jac_mat, 2, max)!=0)==0, 0, 
    sum(apply(jac_mat, 2, max))/sum(apply(jac_mat, 2, max)!=0))
  f <- ifelse(rel * rec == 0, 0, 2 * rel * rec / (rel + rec))
  return(list("rec" = rep(rec, 2), "rel" = rep(rel, 2), "f_score" = rep(f, 2)))
}

# calculate sil score for one view
sil_score <- function(Xinput, row_clustering, col_clustering, index = 2){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    n_views <- length(Xinput)
    n_clusts <- ncol(row_clustering)
    sil_score <- rep(0, length = n_clusts)
    if (index == 2) {
        for (k in 1:n_clusts){
          #select data from specific column clustering
          new_data <- Xinput[, (col_clustering[,k] == 1)]
          spear_dists <- as.matrix(dist(new_data, upper = TRUE, diag = TRUE))
          indices <- row_clustering[, k] == 1
          if ((sum(indices) == 0) |(sum(col_clustering[,k] == 1) == 0) ){
              sil_score[k] <- 0
            }else {
              a_vals <- rowMeans(spear_dists[indices, indices])
              b_vals <- rowMeans(spear_dists[indices, !indices])
              s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
              sil_score[k] <- mean(s_con)
            }
        }
    }
    sil <- ifelse(n_clusts == 1, mean(sil_score),
                 mean(sil_score) - sd(sil_score))
  #return relationships and mean of sil_score across views
  return(list("sil" = sil))
}
sil_score <- function(Xinput, row_clustering, col_clustering, index=2){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    n_clusts <- ncol(row_clustering)
    sil_score <- rep(0, length = n_clusts)
    if (index == 2) {
        clust_one <- col_clustering
        clust_two <- row_clustering
    }else{
        clust_one <- row_clustering
        clust_two <- col_clustering
    }
    for (k in 1:n_clusts){
      #select data from specific column clustering
      if (index == 2){
        new_data <- Xinput[, (clust_one[, k] == 1)]
      }else{
        new_data <- t(Xinput[(clust_one[, k] == 1), ])
      }
      spear_dists <- as.matrix(dist(new_data, upper = TRUE, diag = TRUE))
      indices <- clust_two[, k] == 1
      if ((sum(indices) == 0) || (sum(clust_one[,k] == 1) == 0)){
        sil_score[k] <- 0
      }else {
        a_vals <- rowMeans(spear_dists[indices, indices])
        b_vals <- rowMeans(spear_dists[indices, !indices])
        s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
        sil_score[k] <- mean(s_con)
      }
    }
  if (sum(sil_score) == 0){
    sil <- 0
  }else{
    sil <- ifelse(sum(sil_score != 0) <= 1, sum(sil_score),
                 sum(sil_score) / (sum(sil_score != 0)) -
                 2 * sd(sil_score[sil_score != 0]))
  }
  
  #return relationships and mean of sil_score across views
  return(list("sil" = sil))
}




csr <- function(row_clustering, col_clustering, true_row_clustering, true_col_clustering){
  n_row_cl <- sum(colSums(row_clustering) > 0)
  n_col_cl <- sum(colSums(col_clustering) > 0)
  true_row_cl <- sum(colSums(true_row_clustering) > 0)
  true_col_cl <- sum(colSums(true_col_clustering) > 0)

  #calculate csr
  row_csr <- 1 - abs(n_row_cl - true_row_cl) / (n_row_cl + true_row_cl + 1)
  col_csr <- 1 - abs(n_col_cl - true_col_cl) / (n_col_cl + true_col_cl + 1)
  return(c(row_csr, col_csr))
}

evaluate_simulation_comp <- function(row_clustering, col_clustering, true_row_clustering, true_col_clustering, data_views, index = 2){
  #' row_clust/col_clust: list of the clustering for each view of the rows/columns given by a method
  #' true_row/col_clustering: list of the true clustering for each view of the rows/columns 
  n_views <- length(row_clustering)
  #accuracy table for each view - 1st row is row-clustering, second is column
  csr <- matrix(0, nrow = 2, ncol = n_views)
  rownames(csr) <- c("Row-clustering", "Column-clustering")
  #relevance/recovery for each view
  rel <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rel) <- c("Row-clustering", "Column-clustering")
  rec <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rec) <- c("Row-clustering", "Column-clustering")
  f_score <- matrix(0, nrow = 2, ncol = n_views)
  rownames(f_score) <- c("Row-clustering", "Column-clustering")
  sil_mat <- matrix(0, nrow = 2, ncol = n_views)
  rownames(sil_mat) <- c("Row-clustering", "Column-clustering")
  #go through each view
  for (i in 1:n_views){
    #clust_check <- table(row_clustering[[i]], true_row_clustering[[i]])
    sil_mat[, i] <- rep(sil_score(data_views[[i]], row_clustering[[i]],
                                  col_clustering[[i]], index)$sil, 2)
                    
    jaccard <- jaccard_results(row_clustering[[i]],
                                  col_clustering[[i]],
                                   true_row_clustering[[i]],
                                    true_col_clustering[[i]])
    rel[, i] <- jaccard$rel
    rec[, i] <- jaccard$rec
    f_score[, i] <- jaccard$f_score
    
    #in the results
    csr[, i] <- csr(row_clustering[[i]],
                                  col_clustering[[i]],
                                   true_row_clustering[[i]],
                                    true_col_clustering[[i]])
    
  }

  return(list("CSR" = csr,
              "BiS" = sil_mat,
              "Rel" =  rel,
              "Rec" = rec,
              "F_score" = f_score))
}
