# General solution
update_F <- function(Xinput, Finput, Sinput, Ginput, phi, k){
  #' X: Input matrix
  #' F: row-clustering -- Entire list as input of length n_v
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' k: which view to update
  #' Output: An update for Finput[[k]] 
  
  # Find numerator
  currentF <- Finput[[k]]
  numerator_matrix <- Xinput %*% Ginput %*% t(Sinput)
  denominator_matrix <- currentF %*% Sinput %*% t(Ginput) %*% Ginput %*% t(Sinput)
  #denominator_matrix <- currentF %*% t(currentF) %*% Xinput %*% Ginput %*% t(Sinput)
  numerator <- matrix(0, nrow = nrow(currentF), ncol = ncol(currentF))
  denominator <- matrix(1, nrow = nrow(currentF), ncol = ncol(currentF))
  outputF <- currentF
  temp_den <- 0
  temp_num <- 0
  for (i in 1:nrow(currentF)){
    for (j in 1:ncol(currentF)){
      temp_1 <- 0
      temp_2 <- 0
      temp_3 <- 0
      for (u in 1:length(Finput)){
        if (phi[u,k] !=0){
          temp_1 <- temp_1 + phi[u,k]*Finput[[u]]
          temp_3 <- temp_3 + phi[u,k]
        }
        if (phi[k,u] !=0){
          temp_2 <- temp_2 + phi[k,u]*Finput[[u]]
          temp_3 <- temp_3 + phi[k,u]
        }
      }
      temp_den <- temp_3*Finput[[k]]
      temp_num <- temp_1 + temp_2
      if (is.na(sum(temp_num))){
        print("NA is present")
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else if (sum(temp_num)== 0){
        print("zero sum")
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else{
        numerator[i,j] <- numerator_matrix[i,j] + temp_num[i,j] 
        #numerator is not na
        denominator[i,j] <- denominator_matrix[i,j] + temp_den[i,j]
      }
      outputF[i,j] <- currentF[i,j] * (numerator[i,j]/denominator[i,j])
      if(is.na(outputF[i,j])){print(c(currentF[i,j] , numerator[i,j],denominator[i,j]))}
    }
  }
  return(outputF)
}
update_S <- function(Xinput, Finput, Sinput, Ginput, xi, k){
  #' X: Input matrix
  #' F: row-clustering
  #' S: connection between row-clustering and column-clustering -- Entire list as input of length n_v
  #' G: column-clustering
  #' xi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' k: which view to update
  #' Output: An update for Finput[[k]] 
  
  # Find numerator
  currentS <- Sinput[[k]]
  numerator_matrix <- t(Finput) %*% Xinput %*% Ginput
  denominator_matrix <- t(Finput) %*% Finput %*% currentS %*% t(Ginput) %*% Ginput
  numerator <- matrix(0, nrow = nrow(currentS), ncol = ncol(currentS))
  denominator <- matrix(1, nrow = nrow(currentS), ncol = ncol(currentS))
  outputS <- currentS
  
  for (i in 1:nrow(currentS)){
    for (j in 1:ncol(currentS)){
      temp_1 <- 0
      temp_2 <- 0
      temp_3 <- 0
      for (u in 1:length(Sinput)){
        if (xi[u,k] !=0){
          temp_1 <- temp_1 + xi[u,k]*Sinput[[u]]
          temp_3 <- temp_3 + xi[u,k]
        }
        if (xi[k,u] !=0){
          temp_2 <- temp_2 + xi[k,u]*Sinput[[u]]
          temp_3 <- temp_3 + xi[k,u]
        }
      }
      temp_den <- temp_3*Sinput[[k]]
      temp_num <- temp_1 + temp_2
      if (is.na(sum(temp_num))){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else if (sum(temp_num) == 0){
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else{
        numerator[i,j] <- numerator_matrix[i,j] + temp_num[i,j] 
        denominator[i,j] <- denominator_matrix[i,j] + temp_den[i,j]
      }
      outputS[i,j] <- currentS[i,j] * (numerator[i,j]/denominator[i,j])
    }
  }
  return(outputS)
}
update_G <- function(Xinput, Finput, Sinput, Ginput, psi, k){
  #' X: Input matrix
  #' F: row-clustering
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering -- Entire list as input of length n_v
  #' psi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' k: which view to update
  #' Output: An update for Finput[[k]] 
  
  # Find numerator
  currentG <- Ginput[[k]]
  numerator_matrix <- t(Xinput) %*% Finput %*% Sinput
  denominator_matrix <- currentG %*% t(Sinput) %*% t(Finput) %*% Finput %*% Sinput
  #denominator_matrix <- currentG %*% t(currentG) %*% t(Xinput) %*% Finput %*% Sinput
  numerator <- matrix(0, nrow = nrow(currentG), ncol = ncol(currentG))
  denominator <- matrix(1, nrow = nrow(currentG), ncol = ncol(currentG))
  outputG <- currentG
  temp_den <- 0
  temp_num <- 0
  for (i in 1:nrow(currentG)){
    for (j in 1:ncol(currentG)){
      temp_1 <- 0
      temp_2 <- 0
      temp_3 <- 0
      for (u in 1:length(Ginput)){
        if (psi[u,k] !=0){
          temp_1 <- temp_1 + psi[u,k]*Ginput[[u]]
          temp_3 <- temp_3 + psi[u,k]
        }
        if (psi[k,u] !=0){
          temp_2 <- temp_2 + psi[k,u]*Ginput[[u]]
          temp_3 <- temp_3 + psi[k,u]
        }
      }
      temp_den <- temp_3*Ginput[[k]]
      temp_num <- temp_1 + temp_2
      if (is.na(sum(temp_num))){
        print("na")
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else if (sum(temp_num) == 0){
        #print("zero")
        numerator[i,j] <- numerator_matrix[i,j]  
        denominator[i,j] <- denominator_matrix[i,j]
      } else{
        numerator[i,j] <- numerator_matrix[i,j] + temp_num[i,j] 
        denominator[i,j] <- denominator_matrix[i,j] + temp_den[i,j]
      }
      outputG[i,j] <- currentG[i,j] * (numerator[i,j]/denominator[i,j])
    }
  }
  return(outputG)
}


single_alt_l1_normalisation <- function(Xmatrix){
  #newMatrix <- matrix(0, nrow = nrow(Xmatrix), ncol = ncol(Xmatrix))
  # Introduce Q matrix
  # diagonal matrix of column sums
  # normalises a matrix so that the l1 norm of each column is 1
  Q <- diag(colSums(Xmatrix))
  #solve Q
  newMatrix <- Xmatrix %*% solve(Q)
  return(list("Q" = Q, "newMatrix" = newMatrix))
}

restMultiNMTF_algo <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
  #' Follow this structure:
  #' 1. Normalise X^v s.t ||X||_1 = 1
  #' 2. Initialise F, S, G
  #' 3. Repeat for each view u
  #'    a. Fix ALL, update F^u
  #'    b. Normalise F^u and S^u
  #'    c. Fix ALL, update G^u
  #'    c. Normalise G^u and S^u
  #'    d. Fix ALL, update S^u
  #' until convergence
  #' 
  #' X: Input matrix
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v)
  #' psi: weight on restrictions for S -> matrix of size (n_v x n_v)
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v)
  #' Finit: Inital F matrix
  #' Sinit: Inital S matrix
  #' Ginit: Inital G matrix
  #' Output: Foutput, Soutput, Goutput
  
  # Update view-by-view
  n_v <- length(Xinput)
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  for (v in 1:n_v){
    # Update F
    currentF[[v]] <- update_F(Xinput = Xinput[[v]],
                              Finput = currentF,
                              Sinput = currentS[[v]],
                              Ginput = currentG[[v]],
                              phi = phi, k = v)

    # Normalise F and S
    Q_F <- single_alt_l1_normalisation(currentF[[v]])$Q
    currentF[[v]] <- currentF[[v]] %*% solve(Q_F)
    #if(v==1){print(colSums(currentF[[v]]))}
    #currentF[[v]] <- currentF[[v]] %*% solve(Q_F)
    currentS[[v]] <- Q_F %*% currentS[[v]]

    #if(v==1){
      #print(paste("max relative difference between updates is",max(abs(update_old-update_new))))
      #Q_new <- single_alt_l1_normalisation(update_new)$Q
      #normal_new <- update_new %*% solve(Q_new)
      #print(paste("max relative difference between normalised updates is",max(abs(currentF[[v]]-normal_new))))
    #}
    


    #should be able to replace these lines with 
    #F_normal <- single_alt_l1_normalisation(currentF[[v]])
    #currentF[[v]] <- currentF[[v]] %*% (F_normal$newMatrix)
    #currentS[[v]] <- (F_normal$newMatrix) %*% currentS[[v]]

    # Update G
    currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              psi = psi, k = v)
  
    # Normalise F and S
    Q_G <- single_alt_l1_normalisation(currentG[[v]])$Q
    currentG[[v]] <- currentG[[v]] %*% solve(Q_G)
    # is this transpose needed? surely it's a sqaure matrix
    currentS[[v]] <- currentS[[v]] %*% t(Q_G)
    # Update S
    currentS[[v]] <- update_S(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS,
                              Ginput = currentG[[v]],
                              xi = xi, k = v)
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput,"Soutput" = Soutput,"Goutput" = Goutput))
}

# Main function
restMultiNMTF_run <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
  #' Run Restrictive-Multi NMTF, following the above algorithm
  n_v <- length(Xinput)
  # Normalise Xinput 
  for (v in 1:n_v){
    Xinput[[v]] <- single_alt_l1_normalisation(Xinput[[v]])$newMatrix
  }
  
  # Take Finit, Sinit, Ginit as the initialised latent representations
  currentF <- Finput
  currentS <- Sinput
  currentG <- Ginput
  # Initialising additional parameters
  Xhat <- vector("list", length = n_v)
  total_err <- c()
  # Update until convergence, or for nIter times
  if (is.null(nIter)){
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while (err_diff > 1.0e-4){
      err <- numeric(length = n_v)
      new_parameters <- restMultiNMTF_algo(X = Xinput, 
                                           Finput = currentF, 
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           phi = phi,
                                           xi = xi, 
                                           psi = psi)
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
      #      err_diff <- abs(mean_err - err_temp)
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err-err_temp)
      err_temp <- mean_err
    }
  } else {
    for (t in 1:nIter){
      err <- numeric(length = length(currentF))
      print(t)
      new_parameters <- restMultiNMTF_algo(X = Xinput, 
                                           Finput = currentF, 
                                           Sinput = currentS,
                                           Ginput = currentG,
                                           phi = phi,
                                           xi = xi, 
                                           psi = psi)  
      currentF <- new_parameters$Foutput
      currentS <- new_parameters$Soutput
      currentG <- new_parameters$Goutput
      for (v in 1:n_v){
        Xhat[[v]] <- currentF[[v]] %*% currentS[[v]] %*% t(currentG[[v]])
        err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = total_err))
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
    #print(i)
    #print(row_table[[i]])
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
    #print(column_table[[i]])
    accuracy[1,i] <- accuracyTable(row_table[[i]])
    accuracy[2,i] <- accuracyTable(column_table[[i]])
    adjRandValue[1,i] <- adjustedRandIndex(row_clustering[[i]], true_row_clustering[[i]])
    adjRandValue[2,i] <- adjustedRandIndex(column_clustering[[i]], true_col_clustering[[i]])
    nmiValue[1,i] <- NMI(row_clustering[[i]], true_row_clustering[[i]])
    nmiValue[2,i] <- NMI(column_clustering[[i]], true_col_clustering[[i]])
  }
  return(list("accuracy" = accuracy, "ARI" = adjRandValue, "NMI" = nmiValue))
}


library(clue)
library(mclust)
library(aricode)
## Additional functions
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
