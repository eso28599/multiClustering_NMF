library(mclust)
library(clue)
library(aricode)
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


library(proxy)
internal_met <- function(row_c, col_c, data){
    n_r <- length(unique(row_c))
    n_c <- length(unique(col_c))
    #initialise storage of row/col clusters similarities/distances
    dist_mat_r <- matrix(0, nrow = n_r, ncol = n_c)
    sim_mat_r <- matrix(0, nrow = n_r, ncol = n_c)
    dist_mat_c <- matrix(0, nrow = n_c, ncol = n_r)
    sim_mat_c <- matrix(0, nrow = n_c, ncol = n_r)
    for (i in 1:n_r){
            for (j in 1:n_c){
                sub_mat <- data[row_c == (unique(row_c)[i]), col_c == (unique(col_c)[i])]
                dist_mat_r[i, j] <-  mean(dist(sub_mat))
                sim_mat_r[i, j] <-  mean(simil(sub_mat, method = "cosine"))
                dist_mat_c[j,i] <-  mean(dist(sub_mat, by_rows = FALSE))
                sim_mat_c[j, i] <-  mean(simil(sub_mat, method = "cosine", by_rows = FALSE))
            }
    }
    #find average row
    dist_r <- mean(apply(dist_mat_r, 1, min))# min across a row
    sim_r <- mean(apply(sim_mat_r, 1, max))
    dist_c <- mean(apply(dist_mat_r, 1, min))# min down a col
    sim_c <- mean(apply(sim_mat_r, 1, max))

    return(list("dist" = c(dist_r, dist_c), "sim" = c(sim_r, sim_c)))
}


#just takes the clustering returned by a method
evaluate_simulation <- function(row_clustering, col_clustering, true_row_clustering, true_col_clustering,data_views, error_file){
  #' row_clust/col_clust: list of the clustering for each view of the rows/columns given by a method
  #' true_row/col_clustering: list of the true clustering for each view of the rows/columns 
  #
  n_views <- length(row_clustering)
  
  #list length of number of views
  column_table <- vector("list", length = n_views)
  row_table <- vector("list", length = n_views)
  #accuracy table for each view - 1st row is row-clustering, second is column
  accuracy <- matrix(0, nrow = 2, ncol = n_views)
  rownames(accuracy) <- c("Row-clustering", "Column-clustering")
  #adjusted Rand value for each view
  adjRandValue <- matrix(0, nrow = 2, ncol = n_views)
  rownames(adjRandValue) <- c("Row-clustering", "Column-clustering")
  nmiValue <- matrix(0, nrow = 2, ncol = n_views)
  rownames(nmiValue) <- c("Row-clustering", "Column-clustering")
  #relevance/recovery for each view
  relrev<- matrix(0, nrow = 2, ncol = n_views)
  rownames(adjRandValue) <- c("Row-clustering", "Column-clustering")
  sim_met <- matrix(0, nrow = 2, ncol = n_views)
  rownames(sim_met) <- c("Row-clustering", "Column-clustering")
  dist_met <- matrix(0, nrow = 2, ncol = n_views)
  rownames(dist_met) <- c("Row-clustering", "Column-clustering")
  error_mat <- matrix(0, nrow = 2, ncol = n_views)
  rownames(error_mat) <- c("Row-clustering", "Column-clustering")

  #go through each view
  for (i in 1:n_views){
    #for each view set up a table with the clustering from that view and the true row clustering 
    row_table[[i]] <- table(row_clustering[[i]], true_row_clustering[[i]])
    if (nrow(row_table[[i]]) < ncol(row_table[[i]])){
      #if fewer clusters predicted than in truth, add zeros for those true clusterings
      for (g in 1:(ncol(row_table[[i]]) - nrow(row_table[[i]]))){
        row_table[[i]] <- rbind(row_table[[i]], rep(0,ncol(row_table[[i]])))
      } #if more clusters predicted than in truth, add zeros for those false clusterings
    }else if (nrow(row_table[[i]]) > ncol(row_table[[i]])){
      for (g in 1:(nrow(row_table[[i]]) - ncol(row_table[[i]]))){
        row_table[[i]] <- cbind(row_table[[i]], rep(0,nrow(row_table[[i]])))
      }
    }
    #re order so they match
    row_table[[i]] <- fixTablingClustering(row_table[[i]])
    #check same thing for columns - i.e. same number of predicted and true column cluster
    column_table[[i]] <- table(col_clustering[[i]], true_col_clustering[[i]])
    if (nrow(column_table[[i]]) < ncol(column_table[[i]])){
      for (g in 1:(ncol(column_table[[i]]) - nrow(column_table[[i]]))){
        column_table[[i]] <- rbind(column_table[[i]], rep(0,ncol(column_table[[i]])))
      }
    }else if (nrow(column_table[[i]]) > ncol(column_table[[i]])){
      for (g in 1:(nrow(column_table[[i]]) - ncol(column_table[[i]]))){
        column_table[[i]] <- cbind(column_table[[i]], rep(0,nrow(column_table[[i]])))
      }
    }
    #add fix column clusters
    column_table[[i]] <- fixTablingClustering(column_table[[i]])

    # calculate  results
    accuracy[1, i] <- accuracyTable(row_table[[i]])
    accuracy[2, i] <- accuracyTable(column_table[[i]])
    adjRandValue[1, i] <- adjustedRandIndex(row_clustering[[i]],
                                    true_row_clustering[[i]])
    adjRandValue[2, i] <- adjustedRandIndex(col_clustering[[i]],
                                     true_col_clustering[[i]])
    nmiValue[1, i] <- NMI(as.vector(row_clustering[[i]]),
                                    as.vector(true_row_clustering[[i]]))
    nmiValue[2, i] <- NMI(as.vector(col_clustering[[i]]),
                                     as.vector(true_col_clustering[[i]]))
    #error_mat[, i] <- mean(error_file[[i]][(nIter - 5):nIter])
    error_mat[, i] <- mean(error_file[[1]][(nIter - 10):nIter])
    dist_sim <- internal_met(row_clustering[[i]],
                                  col_clustering[[i]], data_views[[i]])
    dist_met[, i] <- dist_sim$dist
    sim_met[, i] <- dist_sim$sim

  }
  return(list("accuracy" = accuracy,
              "ARI" = adjRandValue, "NMI" = nmiValue, "Sim" = sim_met,
              "Dist" = dist_met, "Error" = error_mat))
}
