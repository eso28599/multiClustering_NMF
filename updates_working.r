# Libraries
library(MASS)
library(eList)

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

#G doesn't work like this because the G's are of different dimensions
#maybe we don't want to put this restriction on the same number of individuals either?
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


a_1<- c(1,2,3)
A_1 <- matrix(c(1,2,3,4),2,2)
# diagonal matrix of column sums
Q_1 <- diag(colSums(A_1))
Q_1_inc <- solve(Q_1)
newMatrix1 <- A_1 %*% Q_1_inc
b_1 <- list(matrix(1,2,2),matrix(1,2,2),matrix(1,2,2))
star_prod(a_1,b_1)


#his method just isn't adding any restrictions if the number of features are the same size
#but then they shouldn't be clustered towards each other anyway!!
#so non-zero psi elements don't make sense if G's aren't the same size
#we can just modify the update step in this case! (same for F and S) make it 
#so that if the matrices are different sizes we don't add them at all
