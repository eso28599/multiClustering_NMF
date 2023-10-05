# Libraries
library(philentropy)
library(MASS)
suppressPackageStartupMessages(library(eList))
library(Matrix)
library(rio)
library(clue)
library(aricode)
library(foreach)
library(parallel)
library(doParallel)
# Get the total number of cores
numberOfCores <- detectCores()
# Register all the cores
registerDoParallel(numberOfCores)

##initialisation
init_mats <- function(X, KK, LL, sigma_I = 0.2){
#' Initialize matrices for multiple data views
#'
#' This function initializes matrices F, S, and G for each data view using singular value decomposition (SVD) and random initialization.
#'
#' @param X A list of matrices representing the data views.
#' @param KK A vector specifying the number of desired columns in the matrix F for each data view.
#' @param LL A vector specifying the number of desired columns in the matrices S and G for each data view.
#' @param sigma_I The standard deviation parameter used in the random initialization of Sinit.
#'
#' @return A list containing the initialized matrices Finit, Ginit, and Sinit.
#'
#' @details
#' For each data view, the function performs the following steps:
#' 1. Performs SVD on the input matrix X[[i]] and extracts the first KK[i] left singular vectors as matrix Finit[[i]].
#' 2. Initializes Sinit[[i]] as the absolute values of the first KK[i] eigenvalues from the diagonal matrix of SVD plus random values sampled from a multivariate normal distribution with mean 0 and covariance matrix sigma_I * identity matrix of size LL[i].
#' 3. Extracts the first LL[i] right singular vectors from the SVD and assigns them to Ginit[[i]].
    
    # Initialisation of F, S, G lists - one matrix for each data view
    n_views <- length(X)
    Finit <- vector("list", length = n_views)
    Sinit <- vector("list", length = n_views)
    Ginit <- vector("list", length = n_views)

    #initialise based on svd
    for(i in 1:n_views) {
        K <- KK[i]
        L <- LL[i]
        ss <- svd(X[[i]])
        Finit[[i]] <- abs(ss$u[, 1:K])
        #using this seems to be better than random - makes sense
        Sinit[[i]] <- abs(diag(ss$d)[1:K, 1:L])
        Sinit[[i]] <- Sinit[[i]] + abs(mvrnorm(n = K,
                   mu = rep(0, K), Sigma = sigma_I * diag(L))[1:K, 1:L])
        Ginit[[i]] <- abs(ss$v[, 1:L])
    }
return(list("Finit" = Finit, "Ginit" = Ginit, "Sinit" = Sinit))
}

## update functions
#define function which takes a list and vector as input 
star_prod <- function(vec,mat_list){ 
  #' mat_list: list of k matrices
  #' vec: a vector of length k
  #' Output: a matrix of same dimensions as the entries in mat_list
  vec_mat <- mclapply(zip(vec, mat_list), function(x) x[[1]] * x[[2]])
  return(Reduce('+', vec_mat))
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
  phi_vec <- (phi + t(phi))[, k]
  if (sum(phi_vec) == 0) {
    outputF <- currentF * ((numerator_matrix) / (denominator_matrix))
  }else {
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
  if (sum(psi) == 0) {
    outputG <- currentG * (numerator_matrix / denominator_matrix)
  }else {
    psi_vec <- (psi + t(psi))[, k]
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
  if (sum(xi) == 0){
    outputS <- currentS * (numerator_matrix  / denominator_matrix)
  }else {
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
  return(list("Foutput" = Foutput, "Soutput" = Soutput, "Goutput" = Goutput))
}

##soft clustering implemented
check_biclusters <- function(Xinput, Foutput, repeats, index){
  # check whether the clusters returned
  # actually correspond to a bicluster or are just noise
  n_views <- length(Xinput)
  n_clusts <- dim(Foutput[[1]])[2]
  # updated results
  scores <- matrix(0, nrow = n_views, ncol = n_clusts)
  thresholds <- get_thresholds(Xinput, Foutput, repeats, index)
  for (i in 1:n_views){
    x_noise <- thresholds$data[[i]]
    for (k in 1:n_clusts){
      x <- Foutput[[i]][, k]
      scores[i, k] <- mean(apply(x_noise, 2,
           function(y) suppressMessages(JSD(rbind(x, y), unit = "log2"))))
    }
  }
  return(list("score" = scores,
          "avg_threshold" = thresholds$avg_score,
          "max_threshold" = thresholds$max_score))
}

get_thresholds <- function(Xinput, Foutput, repeats, index){
    #repeats: minimum value of 2
    n_views <- length(Xinput)
    n_clusts <- dim(Foutput[[1]])[2]
    k_input <- n_clusts * rep(1, length = n_views)
    x_mess <- vector(mode = "list", length = repeats)
    #change into a lapply !!
    for (n in 1:repeats) {
      data_messed <- vector(mode = "list", length = n_views)
      for (i in 1:n_views){
          #correct shuffling
          dims <- dim(Xinput[[i]])
          data_messed[[i]] <- matrix(sample(Xinput[[i]]),
                                  dims[1], dims[2])
      }
      # updated results
      x_mess[[n]] <- restMultiNMTF_run(Xinput = data_messed,
            phi = phi_mat, KK = k_input,  LL = k_input,
             relErr = "Sil", no_clusts = TRUE, stability = FALSE)$Foutput
    }
    avg_score <- c()
    max_score <- c()
    data <- vector(mode = "list", length = n_views)
    for (i in 1:n_views){
      scores <- c()
        for (j in 1:max(repeats - 1,1)){
          data[[i]] <- cbind(data[[i]], x_mess[[j]][[i]])
          for (k in 1:n_clusts){
            #jth repeat, ith view, kth cluster
            x1 <- ((x_mess[[j]])[[i]])[, k]
            for (l in (j + 1):repeats){
              for (m in 1:n_clusts){
                x2 <- ((x_mess[[l]])[[i]])[, m]
                scores <- c(scores,
                   suppressMessages(JSD(rbind(x1, x2), unit = "log2")))
              }
            }
          }
      }
      data[[i]] <- cbind(data[[i]], x_mess[[repeats]][[i]])
      avg_score <- c(avg_score, mean(scores))
      max_score <- c(max_score, max(scores))
    }
  return(list("avg_score" = avg_score,
             "max_score" = max_score, "data" = data))
}


new_sil_score <- function(Xinput, Foutput, row_clustering, col_clustering, index){
    #simultaneously calculates silhouette score for a 
    #clustering as well matching clusters correctly.
    n_views <- length(Xinput)
    n_clusts <- dim(Foutput[[1]])[2]
    relations <- matrix(0, nrow = n_views, n_clusts)
    sil_score <- matrix(0, nrow = n_views, n_clusts)
    for (i in 1:n_views){
      if (index == 2) {
        s_mat <- matrix(0, nrow = n_clusts, n_clusts)
        for (k in 1:n_clusts){
          #select data from specific column clustering
          new_data <- Xinput[[i]][, col_clustering[[i]][,k] == 1]
          spear_dists <- as.matrix(dist(new_data, upper = TRUE, diag = TRUE))
          for (j in 1:n_clusts){
            indices <- row_clustering[[i]][, j] == 1
            a_vals <- rowMeans(spear_dists[indices, indices])
            b_vals <- rowMeans(spear_dists[indices, !indices])
            s_con <- (b_vals - a_vals) / apply(rbind(b_vals, a_vals), 2, max)
            s_mat[k, j] <- mean(s_con)
          }
        }
        relations[i, ] <- apply(s_mat, 1, which.max)
        sil_score[i, ] <- apply(s_mat, 1, max)
      }
      #find the contribution to the s_score of the biclusters from this view
      sil <- mean(sil_score)
    }
  #return relationships and mean of sil_score across views
  return(list("sil" = sil, "scores" = sil_score, "relations" = relations))
}

clustering_res_NMTF <- function(Xinput, Foutput,
               Goutput, Soutput, repeats, index){
  n_views <- length(Foutput)
  row_clustering <- vector("list", length = n_views)
  col_clustering <- vector("list", length = n_views)
  #check clustering and remove if necessary
  biclusts <- check_biclusters(Xinput, Foutput, repeats, index)
  for (i in 1:n_views) {
    # at the minute only matching one - one
    # change to comparing differences etc
    #check if we are matching rows to cols
    row_clustering[[i]] <- apply(Foutput[[i]],
           2, function(x) as.numeric(x > 1 / dim(Foutput[[i]])[1]))
    col_clustering[[i]] <- apply(Goutput[[i]],
           2, function(x) as.numeric(x > 1 / dim(Goutput[[i]])[1]))
  }
  sil_score <- new_sil_score(Xinput, Foutput,
       row_clustering, col_clustering, index)
  s_score <- sil_score$score
  #update realtions and
  #set biclusters that aren't strong enough to 0
  for (i in 1:n_views){
     indices <- ((biclusts$score[i, ]) < biclusts$avg_threshold[i])
     s_score[i, indices] <- 0
     relations <- sil_score$relations[i, ]
     if (index == 1) {
      row_clustering[[i]] <- row_clustering[[i]][, relations]
      row_clustering[[i]][, indices[relations]] <- 0
      col_clustering[[i]][, indices] <- 0
    } else {
      col_clustering[[i]] <- col_clustering[[i]][, relations]
      row_clustering[[i]][, indices] <- 0
      col_clustering[[i]][, indices[relations]] <- 0
    }
  }
  #final sil score
  sil <- sum(s_score) / (sum(s_score != 0))
  #match clusters
  return(list("row_clustering" = row_clustering,
      "col_clustering" = col_clustering, "sil" = sil))
}

restMultiNMTF_main <- function(Xinput, Finput = NULL, Sinput = NULL,
         Ginput = NULL, KK = NULL,
         LL = NULL, phi = NULL, xi = NULL, psi = NULL,
          nIter = 1000, relErr = "Sil", 
         repeats = 5, index = 2, no_clusts = FALSE){
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
  n_v <- length(Xinput)
  # check in correct form
  if (!typeof(Xinput[[1]]) == "double") {
    Xinput <- lapply(Xinput, function(x) as.matrix(x))
    }
  # Normalise Xinput
  Xinput <- lapply(Xinput, function(x) single_alt_l1_normalisation(x)$newMatrix)
  # initialise F, S and G based on svd decomposition if not given
  if (is.null(Finput) | is.null(Ginput) | is.null(Sinput) ) {
        inits <- init_mats(Xinput, KK, LL)
        currentF <- inits$Finit
        currentS <- inits$Sinit
        currentG <- inits$Ginit
  }else {
    # Take Finit, Sinit, Ginit as the initialised latent representations
    currentF <- Finput
    currentS <- Sinput
    currentG <- Ginput
  }
  # Initialising additional parameters
  Xhat <- vector("list", length = n_v)
  # Update until convergence, or for nIter times
  if (is.null(nIter)){
    total_err <- c()
    # Run while-loop until convergence
    err_diff <- 1
    err_temp <- 0
    while ((err_diff > 1.0e-4) & (length(total_err) < 1500)) {
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
        if (relErr == "RelErr") {
          err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]])) / sum(abs(Xinput[[v]]))
        }else {
          err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
        }
      }
      mean_err <- mean(err)
      total_err <- c(total_err, mean_err)
      err_diff <- mean(abs(mean_err - err_temp))
      err_temp <- tail(total_err, n = 3)
    }
  } else {
    total_err <- numeric(length = nIter)
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
        if (relErr == "RelErr") {
          err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]])) / sum(abs(Xinput[[v]]))
        }else{
          err[v] <- sum(abs(Xinput[[v]] - Xhat[[v]]))
        }
      }
      total_err[t] <- mean(err)
    }
  }
  Foutput <- currentF
  Soutput <- currentS
  Goutput <- currentG
  # if only need to obtain factorisation, return values now
  if (no_clusts) {
    return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput))
  }
  print(paste0("got here for", KK[1]))
  # find clustering results and silhouette score
  clusters <- clustering_res_NMTF(Xinput, Foutput,
               Goutput, Soutput, repeats, index)
  print(paste0("got here too for", KK[1]))
  if (is.null(nIter)) {
    error <- mean(tail(total_err, n = 10))
  }else {
    error <- tail(total_err, n = 1)
  }
  return(list("Foutput" = Foutput, "Soutput" = Soutput,
              "Goutput" = Goutput, "Error" = error,
              "All_Error" = total_err, "Sil_score" = clusters$sil,
              "row_clusters" = clusters$row_clustering,
              "col_clusters" = clusters$col_clustering))
}

jaccard <- function(a, b) {
    intersection <- length(intersect(a, b))
    union <-  length(a) + length(b) - intersection
    return(intersection / union)
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
  m <- ncol(row_c) #no of clusters
  n <- ncol(true_r)
  samps <- 1:nrow(row_c)
  feats <- 1:nrow(col_c)
  #initialise storage of jaccard index between pairs
  jac_mat <- matrix(nrow = m, ncol = n)
  for (i in 1:m){
    r_i <- samps[row_c[, i] == 1]
    c_i <- feats[col_c[, i] == 1]
    m_i <- cart_prod(r_i, c_i)
    for (j in 1:n){
        tr_i <- samps[true_r[, j] == 1]
        tc_i <- samps[true_c[, j] == 1]
        m_j <- cart_prod(tr_i, tc_i)
        jac_mat[i, j] <- jaccard(m_i, m_j)
    }
  }
  if (stability) {
    return(apply(jac_mat, 2, max))
  }
  rel <- mean(apply(jac_mat, 1, max))
  rev <- mean(apply(jac_mat, 2, max))
  f <- 2 * rel * rev / (rel + rev)
  return(list("rev" = rep(rev, 2), "rel" = rep(rel, 2), "f_score" = rep(f, 2)))
}

stability_check <- function(Xinput, results,
                     k, phi, xi, psi, nIter, relErr,
                     repeats, index, no_clusts, sample_rate = 0.9,
                     n_stability = 5, stab_thres = 0.2){

    n_views <- length(Xinput)
    # initialise storage of results
    jacc <- matrix(0, nrow = n_views, ncol = k[1])
    jacc_rand <- matrix(0, nrow = n_views, ncol = k[1])
    for (t in 1:n_stability){
      new_data <- vector(mode = "list", length = n_views)
      row_samples <- vector(mode = "list", length = n_views)
      col_samples <- vector(mode = "list", length = n_views)
      #turn this into a function to be used with lapply
      for(i in 1:n_views){
        dims <- dim(Xinput[[i]])
        row_samples[[i]] <- sample(dims[1], (dims[1] * sample_rate))
        col_samples[[i]] <- sample(dims[2], (dims[2] * sample_rate))
        new_data[[i]] <- Xinput[[i]][row_samples[[i]], col_samples[[i]]]
      }

      new_results <- restMultiNMTF_main(new_data, Finput = NULL, Sinput = NULL,
          Ginput = NULL, k,
          k, phi, xi, psi,
            nIter, relErr = "Sil",
          repeats, index = 2)
      #compare results
      #extract results
      for(i in 1:n_views){
        jacc[i, ] <- jacc[i, ] + jaccard_results(new_results$row_clusters[[i]],
                 new_results$col_clusters[[i]],
              results$row_clusters[[i]][row_samples[[i]], ],
               results$col_clusters[[i]][col_samples[[i]], ], TRUE)
        jacc_rand[i, ] <- jacc_rand[i, ] + jaccard_results(apply(new_results$row_clusters[[i]], 2, sample),
                 apply(new_results$col_clusters[[i]], 2, sample),
              results$row_clusters[[i]][row_samples[[i]], ],
               results$col_clusters[[i]][col_samples[[i]],], TRUE)
      }
    }
    jacc <- jacc / n_stability
    jacc_rand <- jacc_rand / n_stability
    for (i in 1:n_views){
      #set clusters not deemed stable to have 0 members
      results$row_clusters[[i]][, jacc[i, ] < jacc_rand[i, ]] <- 0
      results$col_clusters[[i]][, jacc[i, ] < jacc_rand[i, ]] <- 0
    }
    return(results)
}
init_rest_mats <- function(mat, n_v){
  if (is.null(mat)){
    return(matrix(0, nrow = n_v, ncol = n_v))
  }else{
    return(0.01 * mat)
  }
}

restMultiNMTF_run <- function(Xinput, Finput=NULL, Sinput=NULL, 
            Ginput=NULL, KK=NULL, LL=NULL, phi=NULL, xi=NULL, psi=NULL, 
            nIter=1000, relErr=TRUE, k_max =6, repeats = 5, index = 2, no_clusts = FALSE, 
             sample_rate = 0.8, n_stability = 5, stab_thres = 0.8, stability = TRUE){

  #' @param k_max integer, default is 6, must be greater than 2, largest value of k to be considered initially,
  # initialise phi etc matrices as zeros if not specified
  # otherwise multiply by given parameter
  n_v <- length(Xinput)
  # initialise restriction matrices if not specified 
  # views with no restrictions require no input
  phi <- init_rest_mats(phi, n_v)
  psi <- init_rest_mats(psi, n_v)

  # if number of clusters has been specified method can be applied straight away
  if ((!is.null(KK)) || (!is.null(LL))) {
    results <- restMultiNMTF_main(Xinput, Finput, Sinput, Ginput,
                     KK, LL, phi, xi, psi, nIter, relErr,
                      repeats, index, no_clusts)
    # if using the original data, we want to perform stability analysis 
    # otherwise we want the results
    if (stability) {
      return(stability_check(Xinput, results,
                     KK, phi, xi, psi, nIter, relErr,
                    repeats, index, no_clusts, sample_rate,
                     n_stability, stab_thres))
    }else {
      return(results)
    }
  }
  # define set of k_s to consider
  KK <- 2:k_max
  k_vec <- rep(1, n_v)
  #initialise storage of results
  res_list <- vector("list", length = length(KK))
  err_list <- rep(0, length = length(KK))

  # apply method for each k to be considered
  foreach(i = 1:length(KK)) %dopar% {
    k <- KK[i] * k_vec
    print(paste0(KK[i], " started"))
    res_list[[i]] <- restMultiNMTF_main(Xinput, Finput, Sinput, Ginput,
                     k, k, phi, xi, psi, nIter, relErr,
                     repeats, index, no_clusts)
    if(relErr == "Sil"){
      err_list[i] <-res_list[[i]][["Sil_score"]]
    }else{
      err_list[i] <- res_list[[i]][["Error"]]
    }
    print(paste0(KK[i], " done"))
  }
  #find value of k of lowest error
  test <- KK[ifelse(relErr == "Sil", which.max(err_list), which.min(err_list))]
  max_i <- 6
  # if best performing k is the largest k considered
  # apply method to k + 1 until this is no longer the case
  while (test == max_i) {
    max_i <- max_i + 1
    KK <- c(KK, max_i)
    k <- max_i * k_vec
    res_list[[max_i]] <- restMultiNMTF_main(Xinput, Finput, Sinput, Ginput,
                     k, k, phi, xi, psi, nIter, relErr,
                     repeats, index, no_clusts)
    if(relErr == "Sil"){
      err_list[max_i] <- res_list[[max_i]][["Sil_score"]]
      test <- KK[which.max(err_list)]
    }else {
      err_list[max_i] <- res_list[[max_i]][["Error"]]
      test <- KK[which.min(err_list)]
    }
  }
  if(relErr == "Sil"){
    k <- which.max(err_list)
    results <- res_list[[k]]
    k_vec <- rep(KK[k], length = n_v)
    return(stability_check(Xinput, results,
                     k_vec, phi, xi, psi, nIter, relErr,
                    repeats, index, no_clusts, sample_rate, n_stability, stab_thres))
  }else{
    k <- which.min(err_list)
    results <- res_list[[k]]
    k_vec <- rep(KK[k], length = n_v)
    return(stability_check(Xinput, results,
                     k_vec, phi, xi, psi, nIter, relErr,
                    repeats, index, no_clusts, sample_rate, n_stability, stab_thres))
  }
}
