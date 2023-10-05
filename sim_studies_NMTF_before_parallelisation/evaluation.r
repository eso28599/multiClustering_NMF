library(mclust)

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

evaluate_simulation <- function(X_nmtf, true_row_clustering, true_col_clustering){
  #' X_nmtf: Output of restMultiNMTF_run
  #' true_row/col_clustering: What the name states
  #' 
  
  #list length of number of views
  row_clustering <- vector("list", length = length(X_nmtf$Foutput))
  column_clustering <- vector("list", length = length(X_nmtf$Foutput))
  for (i in 1:length(X_nmtf$Foutput)){
    row_clustering[[i]] <- apply(X_nmtf$Foutput[[i]], 1, which.max)
    column_clustering[[i]] <- apply(X_nmtf$Goutput[[i]], 1, which.max)
  }
  #list length of number of views
  column_table <- vector("list", length = length(X_nmtf$Foutput))
  row_table <- vector("list", length = length(X_nmtf$Foutput))
  #accuracy table for each view - 1st row is row-clustering, second is column
  accuracy <- matrix(0, nrow = 2, ncol = length(X_nmtf$Foutput))
  rownames(accuracy) <- c("Row-clustering", "Column-clustering")
  #adjusted Rand value for each view
  adjRandValue <- matrix(0, nrow = 2, ncol = length(X_nmtf$Foutput))
  rownames(adjRandValue) <- c("Row-clustering", "Column-clustering")
  nmiValue <- matrix(0, nrow = 2, ncol = length(X_nmtf$Foutput))
  rownames(nmiValue) <- c("Row-clustering", "Column-clustering")

  #go through each view
  for (i in 1:length(X_nmtf$Foutput)){
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
    column_table[[i]] <- table(column_clustering[[i]], true_col_clustering[[i]])
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
    accuracy[1,i] <- accuracyTable(row_table[[i]])
    accuracy[2,i] <- accuracyTable(column_table[[i]])
    adjRandValue[1,i] <- adjustedRandIndex(row_clustering[[i]], true_row_clustering[[i]])
    adjRandValue[2,i] <- adjustedRandIndex(column_clustering[[i]], true_col_clustering[[i]])
    nmiValue[1,i] <- NMI(row_clustering[[i]], true_row_clustering[[i]])
    nmiValue[2,i] <- NMI(column_clustering[[i]], true_col_clustering[[i]])
  }
  return(list("accuracy" = accuracy, "ARI" = adjRandValue, "NMI" = nmiValue))
}




#modification - just takes the clustering returned by a method
evaluate_simulation <- function(row_clustering,col_clustering, true_row_clustering, true_col_clustering){
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
    column_table[[i]] <- table(column_clustering[[i]], true_col_clustering[[i]])
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
    accuracy[1,i] <- accuracyTable(row_table[[i]])
    accuracy[2,i] <- accuracyTable(column_table[[i]])
    adjRandValue[1,i] <- adjustedRandIndex(row_clustering[[i]], true_row_clustering[[i]])
    adjRandValue[2,i] <- adjustedRandIndex(column_clustering[[i]], true_col_clustering[[i]])
    nmiValue[1,i] <- NMI(row_clustering[[i]], true_row_clustering[[i]])
    nmiValue[2,i] <- NMI(column_clustering[[i]], true_col_clustering[[i]])
    
    # we also want 

  }
  return(list("accuracy" = accuracy, "ARI" = adjRandValue, "NMI" = nmiValue))
}




