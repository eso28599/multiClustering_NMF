# Libraries
library(MASS)
library(eList)
library(MASS)
library(Matrix)
library(openxlsx)
library(rio)
library(clue)
library(mclust)
library(aricode)


#define function which takes a list and vector as input 
star_prod <- function(vec,mat_list){ 
  #' mat_list: list of k matrices
  #' vec: a vector of length k
  #' Output: a matrix of same dimensions as the entries in 
  vec_mat <- lapply(zip(vec,mat_list),function(x) x[[1]]*x[[2]])
  return(Reduce('+',vec_mat))
}
#updates F for a specified dataview 
update_F <- function(Xinput, Finput, Sinput, Ginput, phi, k){
  #' X: Input matrix
  #' F: row-clustering -- Entire list as input of length n_v
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering
  #' phi: weight on restrictions for F -> matrix of size (n_v x n_v) - a sparse matrix
  #' k: which view to update
  #' Output: An update for Finput[[k]] 
  
  # Find numerator
  currentF <- Finput[[k]]
  numerator_matrix <- Xinput %*% Ginput %*% t(Sinput)
  denominator_matrix <- currentF %*% Sinput %*% t(Ginput) %*%  Ginput %*% t(Sinput)

  #calculate the column vector based on phi that is needed
  phi_vec <- (phi+t(phi))[,k]
  if(sum(phi_vec)==0){
    outputF <- currentF * ((numerator_matrix) / (denominator_matrix))
  }else{
    num_mat_prod <- star_prod(phi_vec, Finput)
    denom_mat_prod <- sum(phi_vec) * currentF
    outputF <- currentF * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod))
  }
  return(outputF)
}

update_G <- function(Xinput, Finput, Sinput, Ginput, psi, k){
  #' X: Input matrix
  #' F: row-clustering 
  #' S: connection between row-clustering and column-clustering
  #' G: column-clustering -- Entire list as input of length n_v
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v) - a sparse matrix
  #' k: which view to update
  #' Output: An update for Ginput[[k]]
  
  # Find numerator
  currentG <- Ginput[[k]]
  numerator_matrix <- t(Xinput) %*% Finput %*% Sinput
  denominator_matrix <- currentG  %*% t(Sinput) %*% t(Finput) %*% Finput %*% Sinput
  #check if they are all the same size

  #calculate the column vector based on phi that is needed
  if(sum(psi)==0){
    outputG <- currentG * (numerator_matrix / denominator_matrix)
  }else{
    psi_vec <- (psi + t(psi))[,k]
    num_mat_prod <- star_prod(psi_vec, Ginput)
    denom_mat_prod <- sum(psi_vec) * currentG
    outputG <- currentG * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod))
  }
  return(outputG)
}

update_S <- function(Xinput, Finput, Sinput, Ginput, xi, k){
  #' X: Input matrix
  #' F: row-clustering
  #' S: connection between row-clustering and column-clustering -- Entire list as input of length n_v
  #' G: column-clustering
  #' xi: weight on restrictions for G -> matrix of size (n_v x n_v) - a sparse matrix
  #' k: which view to update
  #' Output: An update for Sinput[[k]]
  
  # Find numerator
  currentS <- Sinput[[k]]
  numerator_matrix <- t(Finput) %*% Xinput %*% Ginput
  denominator_matrix <- t(Finput) %*% Finput %*% currentS %*% t(Ginput) %*% Ginput

  #calculate the column vector based on phi that is needed
  if(sum(xi)==0){
    outputS <- currentS * (numerator_matrix  / denominator_matrix)
  }else{
    xi_vec <- (xi + t(xi))[, k]
    num_mat_prod <- star_prod(xi_vec, Sinput)
    denom_mat_prod <- sum(xi_vec) * currentS
    outputS <- currentS * ((numerator_matrix + num_mat_prod) / (denominator_matrix + denom_mat_prod))
  }
  
  return(outputS)
}

single_alt_l1_normalisation <- function(Xmatrix){
  # Introduce Q matrix
  # diagonal matrix of column sums
  # normalises a matrix so that the l1 norm of each column is 1
  Q <- diag(colSums(Xmatrix))
  #solve Q
  newMatrix <- Xmatrix %*% solve(Q)
  return(list("Q" = Q, "newMatrix" = newMatrix))
}

restMultiNMTF_algo <- function(Xinput, Finput, Sinput, Ginput, phi, xi, psi, nIter){
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
    F_normal <- single_alt_l1_normalisation(currentF[[v]])
    currentF[[v]] <- F_normal$newMatrix
    currentS[[v]] <- (F_normal$Q) %*% currentS[[v]]

    # Update G
    currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              psi = psi, k = v)
  
    # Normalise G and S
    G_normal <- single_alt_l1_normalisation(currentG[[v]])
    currentG[[v]] <- G_normal$newMatrix
    currentS[[v]] <- currentS[[v]] %*% G_normal$Q
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
restMultiNMTF_run <- function(Xinput, Finput=NULL, Sinput=NULL, Ginput=NULL, KK=NULL, LL=NULL, phi=NULL, xi=NULL, psi=NULL, nIter=NULL){
  #' Run Restrictive-Multi NMTF, following the above algorithm
  #' 
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
  n_v <- length(Xinput)
  # Normalise Xinput
  Xinput <- lapply(Xinput,function(x) single_alt_l1_normalisation(x)$newMatrix)
  # initialise phi etc matrices as zeros if not specified 
  if(is.null(phi)){ phi <- matrix(0, nrow = n_v, ncol = n_v)}
  if(is.null(psi)){ psi <- matrix(0, nrow = n_v, ncol = n_v)}
  if(is.null(xi)){ xi <- matrix(0, nrow = n_v, ncol = n_v)}

  # initialise F, S and G based on svd decomposition if not given
  if(is.null(Finput)|is.null(Ginput)|is.null(Sinput)){
        inits <- init_mats(Xinput,KK,LL)
        currentF <- inits$Finit
        currentS <- inits$Sinit
        currentG <- inits$Ginit
  }else{
    # Take Finit, Sinit, Ginit as the initialised latent representations
    currentF <- Finput
    currentS <- Sinput
    currentG <- Ginput
  }
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
      total_err <- c(total_err, mean_err)
      err_diff <- abs(mean_err-err_temp)
      err_temp <- mean_err
    }
  } else {
    for (t in 1:nIter){
      err <- numeric(length = length(currentF))
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

#initialisation function
init_mats <- function(X,KK,LL){
    # Initialisation of F, S, G lists - one matrix for each data view
    n_views <- length(X)
    Finit <- vector("list", length = n_views)
    Sinit <- vector("list", length = n_views)
    Ginit <- vector("list", length = n_views)

    
    #initialise based on svd - why?
    for (i in 1:n_views){
        K <- KK[i]
        L <- LL[i]
        ss <- svd(X[[i]])
        Finit[[i]] <- abs(ss$u[,1:K])
        #using this seems to be better than random - makes sense
        #Sinit[[i]] <- abs(diag(ss$d)[1:K,1:L])
        Sinit[[i]] <- abs(mvrnorm(n = K, mu = rep(0,K), Sigma = diag(L))[1:K,1:L])
        Ginit[[i]] <- abs(ss$v[,1:L])
}
return(list("Finit"=Finit , "Ginit"=Ginit , "Sinit"=Sinit))
}

##### synthetic data generation

#preparation of linked dataviews
# function to generate a dataview for specified bicluster dimensions
# same number of row and column clusters 
# row clusters remain common across all views, columns change 

one_view <- function(row_dims,col_dims,noise=10){
    #'row_dims: vector of sizes of each row cluster in this view
    #'col_dims: vector of sizes of each row cluster in this view
    #'noise: variance of the noise added to views
    #' 
    #'view: matrix of this dataview
    #'truth_row: vector indicating membership of row clusters
    #'truth_col: vector indicating membership of col clusters
    n_row <- length(row_dims)
    n_col <- length(col_dims)


    # Create a list as input for block-diagonal data generation
    # list length of no of clusters in this first view
    inputList <- vector("list", length = max(n_row, n_col))
    # generate a mvn with n=number of individuals of view and no of features
    # equal to number of features of the view
    # mean of 5 for each column
    # covariance matrix is identity - each feature is independent 
    #define vectors indicating true row/column cluster membership for each view
    true_row <- c()
    true_col <- c()
    for(i in 1: n_row ){
        inputList[[i]] <- mvrnorm(n =row_dims[i], mu = rep(5,col_dims[i]), Sigma = diag(col_dims[i]))
        true_row <- c(true_row,rep(i,row_dims[i]))
        true_col <- c(true_col,rep(i,col_dims[i]))
    }

    X_view <- as.matrix(bdiag(inputList)) 
    # generate noise, MVN with no correlation, zero mean and variance of 10
    X_noise <- mvrnorm(n = nrow(X_view), mu = rep(0, ncol(X_view)), Sigma = diag(noise,ncol(X_view)))
    # add noise to first view
    X_view <- X_view + X_noise 

    return(list(view= X_view, truth_row = true_row, truth_col = true_col))

}


multi_view <- function(rowClusters,colClusters,noise=10,seed=FALSE){
    #'rowClusters: n length list of vectors of row cluster sizes in each view
    #'rowClusters: n length list of vectors of column cluster sizes in each view
    #' seed: logical indicator, default is FALSE, if true sets seed so same data is generated each time
    #' 
    #'data_views: n length list of data views
    #'truth_rows: n length list of vectors of true row cluster membership for each view
    #'col_rows: n length list of vectors of true column cluster membership for each view
    
    if(seed){
        set.seed(20)
    }
    n_views <- length(rowClusters)
    # Introduce simulated data-views- where we store the views
    X_trial <- vector("list", length = n_views)
    #list to store row clusterings for each dataview
    true_row_clusterings <- vector("list", length = n_views) 
    true_col_clusterings <- vector("list", length = n_views) 

    #generate data and store
    for(i in 1:n_views){
        data_i <- one_view(rowClusters[[i]],colClusters[[i]],noise)
        X_trial[[i]] <- data_i$view
        true_row_clusterings[[i]] <- data_i$truth_row
        true_col_clusterings[[i]] <- data_i$truth_col
    }

    return(list(data_views= X_trial, truth_rows = true_row_clusterings, truth_cols = true_col_clusterings))
}
