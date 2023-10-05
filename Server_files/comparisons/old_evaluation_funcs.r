suppressPackageStartupMessages(library(mclust))
library(clue)
library(aricode)
library(cluster)
suppressPackageStartupMessages(library(proxy))

internal_met <- function(row_c, col_c, data){
    sil_r <- summary(silhouette(row_c, stats::dist(data)))
    sil_c <- summary(silhouette(col_c, stats::dist(t(data))))
    if(length(sil_r) == 2){
      sil_r <- -1
    }else{
      sil_r <- sil_r$avg.width
    }
    if(length(sil_c) == 2){
      sil_c <- -1
    }else{
      sil_c <- sil_c$avg.width
    }
    return(list("sil" = c(sil_r, sil_c)))
}
#performance measures document
fixTablingClustering <- function(tableInput, max = TRUE){
  if (nrow(tableInput) != ncol(tableInput)){
    stop("PLease enter a square matrix")
  }
  lengthTable <- nrow(tableInput)
  orderTemp <- solve_LSAP(tableInput, maximum = max)
  tableTemp <- tableInput
  for (i in 1:lengthTable){
    tableTemp[i,] <- tableInput[which(orderTemp==i),]
  }
  return(tableTemp)
}

# Measuer Accuracy
accuracyTable <- function(tableInput){
  sumDiag <- sum(diag(tableInput))
  sumAll <- sum(tableInput)
  acc <- sumDiag/sumAll
  return(acc)
}

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <-  length(a) + length(b) - intersection
    return(intersection / union)
}

cart_prod <- function(a, b) {
  prod <- c()
  for (k in 1:length(a)){
    prod <- c(prod, paste(a[k], b))
  }
  return(prod)
}
jaccard_results <- function(row_c, col_c, true_r, true_c){
  clusts <- unique(row_c)
  true_clusts <- unique(true_r)
  m <- length(clusts)
  n <- length(true_clusts)
  samps <- 1:length(row_c)
  feats <- 1:length(col_c)
  #initialise jaccard index between pairs
  jac_mat <- matrix(nrow = m, ncol = n)
  for (i in 1:m){
    r_i <- samps[row_c == clusts[i]]
    c_i <- feats[col_c == clusts[i]]
    m_i <- cart_prod(r_i, c_i)
    for (j in 1:n){
        tr_i <- samps[true_r == true_clusts[j]]
        tc_i <- samps[true_c == true_clusts[j]]
        m_j <- cart_prod(tr_i, tc_i)
        jac_mat[i, j] <- jaccard(m_i, m_j)
    }
  }
  print(jac_mat)
  rel <- mean(apply(jac_mat, 1, max))
  rev <- mean(apply(jac_mat, 2, max))
  f <- 2 * rel * rev / (rel + rev)
  return(list("rev" = rep(rev, 2), "rel" = rep(rel, 2), "f_score" = rep(f, 2)))
}

evaluate_simulation_all <- function(row_clustering, col_clustering, true_row_clustering, true_col_clustering,data_views){
  #' row_clust/col_clust: list of the clustering for each view of the rows/columns given by a method
  #' true_row/col_clustering: list of the true clustering for each view of the rows/columns 
  n_views <- length(row_clustering)
  #accuracy table for each view - 1st row is row-clustering, second is column
  accuracy <- matrix(0, nrow = 2, ncol = n_views)
  rownames(accuracy) <- c("Row-clustering", "Column-clustering")
  #adjusted Rand value for each view
  adjRandValue <- matrix(0, nrow = 2, ncol = n_views)
  rownames(adjRandValue) <- c("Row-clustering", "Column-clustering")
  nmiValue <- matrix(0, nrow = 2, ncol = n_views)
  rownames(nmiValue) <- c("Row-clustering", "Column-clustering")
  #relevance/recovery for each view
  rel <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rel) <- c("Row-clustering", "Column-clustering")
  rev <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rev) <- c("Row-clustering", "Column-clustering")
  f_score <- matrix(0, nrow = 2, ncol = n_views)
  rownames(f_score) <- c("Row-clustering", "Column-clustering")
  sil_mat <- matrix(0, nrow = 2, ncol = n_views)
  rownames(sil_mat) <- c("Row-clustering", "Column-clustering")
  #go through each view
  for (i in 1:n_views){
    clust_check <- table(row_clustering[[i]], true_row_clustering[[i]])
    adjRandValue[1, i] <- adjustedRandIndex(row_clustering[[i]],
                                    true_row_clustering[[i]])
    adjRandValue[2, i] <- adjustedRandIndex(col_clustering[[i]],
                                     true_col_clustering[[i]])
    nmiValue[1, i] <- NMI(as.vector(row_clustering[[i]]),
                                    as.vector(true_row_clustering[[i]]))
    nmiValue[2, i] <- NMI(as.vector(col_clustering[[i]]),
                                     as.vector(true_col_clustering[[i]]))
    sil_mat[, i] <- internal_met(row_clustering[[i]],
                                  col_clustering[[i]], data_views[[i]])$sil
    jaccard <- jaccard_results(row_clustering[[i]],
                                  col_clustering[[i]],
                                   true_row_clustering[[i]],
                                    true_col_clustering[[i]])
    rel[, i] <- jaccard$rel
    rev[, i] <- jaccard$rev
    f_score[, i] <- jaccard$f_score
    accuracy[, i] <- ifelse(nrow(clust_check) == ncol(clust_check), 1, 0)
  }

  return(list("Accuracy" = accuracy,
              "ARI" = adjRandValue, "NMI" = nmiValue,
              "Sil" = sil_mat, 
              "Rel" =  rel, 
              "Rev" = rev, 
              "F_score" = f_score))
}

evaluate_simulation_comp <- function(row_clustering, col_clustering, true_row_clustering, true_col_clustering, data_views){
  #' row_clust/col_clust: list of the clustering for each view of the rows/columns given by a method
  #' true_row/col_clustering: list of the true clustering for each view of the rows/columns 
  n_views <- length(row_clustering)
  #accuracy table for each view - 1st row is row-clustering, second is column
  accuracy <- matrix(0, nrow = 2, ncol = n_views)
  rownames(accuracy) <- c("Row-clustering", "Column-clustering")
  #relevance/recovery for each view
  rel <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rel) <- c("Row-clustering", "Column-clustering")
  rev <- matrix(0, nrow = 2, ncol = n_views)
  rownames(rev) <- c("Row-clustering", "Column-clustering")
  f_score <- matrix(0, nrow = 2, ncol = n_views)
  rownames(f_score) <- c("Row-clustering", "Column-clustering")
  sil_mat <- matrix(0, nrow = 2, ncol = n_views)
  rownames(sil_mat) <- c("Row-clustering", "Column-clustering")
  #go through each view
  for (i in 1:n_views){
    clust_check <- table(row_clustering[[i]], true_row_clustering[[i]])
    sil_mat[, i] <- internal_met(row_clustering[[i]],
                                  col_clustering[[i]], data_views[[i]])$sil
    jaccard <- jaccard_results(row_clustering[[i]],
                                  col_clustering[[i]],
                                   true_row_clustering[[i]],
                                    true_col_clustering[[i]])
    rel[, i] <- jaccard$rel
    rev[, i] <- jaccard$rev
    f_score[, i] <- jaccard$f_score
    accuracy[, i] <- ifelse(nrow(clust_check) == ncol(clust_check), 1, 0)
  }

  return(list("Accuracy" = accuracy,
              "Sil" = sil_mat,
              "Rel" =  rel,
              "Rev" = rev,
              "F_score" = f_score))
}
