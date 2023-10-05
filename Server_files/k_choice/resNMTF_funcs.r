# Libraries
library(MASS)
library(eList)
library(Matrix)
library(rio)
library(clue)
library(aricode)


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
        Sinit[[i]] <- abs(diag(ss$d)[1:K,1:L])
        #Sinit[[i]] <- abs(mvrnorm(n = K, mu = rep(0,K), Sigma = diag(L))[1:K,1:L])
        Ginit[[i]] <- abs(ss$v[,1:L])
}
return(list("Finit"=Finit , "Ginit"=Ginit , "Sinit"=Sinit))
}

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

# function which determines clustering by NMTF
#list length of number of views
clustering_res_NMTF <- function(Foutput,Goutput,n_views){
  row_clustering <- vector("list", length = n_views )
  column_clustering <- vector("list", length = n_views)
  for (i in 1:n_views){
    row_clustering[[i]] <- apply(Foutput[[i]], 1, which.max)
    column_clustering[[i]] <- apply(Goutput[[i]], 1, which.max)
  }
  return(list("row_clustering" = row_clustering, "col_clustering" = column_clustering))
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
  # check in correct form
  if (!typeof(Xinput[[1]]) == "double") {
    Xinput <- lapply(Xinput, function(x) as.matrix(x))
    }
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
  total_err <- numeric(length = nIter)
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
      total_err[t] <- mean(err)
      #mean_err <- mean(err)
      #total_err <- c(total_err, mean_err)
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  # also add in clustering results
  clusters <- clustering_res_NMTF(Foutput,Goutput,n_v)

  return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = total_err,
              "row_clusters" = clusters$row_clustering,
              "col_clusters" = clusters$col_clustering))
}
