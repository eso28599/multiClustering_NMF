#removes col of S
clustering_res_NMTF2 <- function(Foutput, Goutput, Soutput, n_views){
  row_clustering <- vector("list", length = n_views )
  column_clustering <- vector("list", length = n_views)
  n_clusts <- dim(Soutput[[1]])[1]
  for (i in 1:n_views){
    #update G
    update_mat <- matrix(data = 0, nrow = n_clusts, ncol = n_clusts)
    for (j in 1:n_clusts){
      match <- which.max(Soutput[[i]][, j])
      Soutput[[i]][match, ] <- 0 # update so not considered next time
      update_mat[match, j] <- 1
    }
    row_clustering[[i]] <- apply(Foutput[[i]], 1, which.max)
    column_clustering[[i]] <- apply(Goutput[[i]] %*% update_mat, 1, which.max)
  }
  return(list("row_clustering" = row_clustering, "col_clustering" = column_clustering))
}

#normalises differently
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
    SG_normal <- single_alt_l1_normalisation(t(currentS[[v]] %*% t(currentG[[v]])))
    #F_normal <- single_alt_l1_normalisation(currentF[[v]])
    #currentF[[v]] <- F_normal$newMatrix
    currentF[[v]] <- currentF[[v]] %*% (SG_normal$Q)
    currentS[[v]] <- solve(SG_normal$Q) %*% currentS[[v]]

    # Update G
    currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              psi = psi, k = v)
  
    # Normalise G and S
    FS_normal <- single_alt_l1_normalisation(currentF[[v]] %*% currentS[[v]])
    #G_normal <- single_alt_l1_normalisation(currentG[[v]])
    currentG[[v]] <- currentG[[v]] %*% (FS_normal$Q)
    currentS[[v]] <- currentS[[v]] %*% solve(FS_normal$Q)
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

#
clustering_res_NMTF_update <- function(Foutput, Goutput, Soutput, n_views, w){
  row_clustering <- vector("list", length = n_views )
  column_clustering <- vector("list", length = n_views)
  n_clusts <- dim(Soutput[[1]])[1]
  cut_off <- (1/n_clusts + w)/2
  for (i in 1:n_views){
    #update G
    update_mat <- matrix(data = 0, nrow = n_clusts, ncol = n_clusts)
    for (j in 1:n_clusts){
      match <- which.max(Soutput[[i]][, j])
      Soutput[[i]][match, ] <- 0 # update so not considered next time
      update_mat[match, j] <- 1
    }
    row_clustering[[i]] <- t(apply(Foutput[[i]]/ apply(nmtf_results$Foutput[[i]],1,max),
           1, function(x) as.numeric(x > cut_off)))
    column_clustering[[i]] <- t(apply((Goutput[[i]] / apply(nmtf_results$Goutput[[i]],1,max)) %*% update_mat,
     1, function(x) as.numeric(x > cut_off)))
  }
  return(list("row_clustering" = row_clustering, "col_clustering" = column_clustering))
}

clustering_res_NMTF <- function(Foutput, Goutput, Soutput, n_views, w){
  row_clustering <- vector("list", length = n_views )
  column_clustering <- vector("list", length = n_views)
  n_clusts <- dim(Soutput[[1]])[1]
  cut_off <- (1/n_clusts + w)/2
  for (i in 1:n_views){
    #update G
    update_mat <- matrix(data = 0, nrow = n_clusts, ncol = n_clusts)
    for (j in 1:n_clusts){
      match <- which.max(Soutput[[i]][, j])
      Soutput[[i]][match, ] <- 0 # update so not considered next time
      update_mat[match, j] <- 1
    }
    row_clustering[[i]] <- t(apply(Foutput[[i]]/ apply(nmtf_results$Foutput[[i]],1,max),
           1, function(x) as.numeric(x > cut_off)))
    column_clustering[[i]] <- t(apply((Goutput[[i]] / apply(nmtf_results$Goutput[[i]],1,max)) %*% update_mat,
     1, function(x) as.numeric(x > cut_off)))
  }
  return(list("row_clustering" = row_clustering, "col_clustering" = column_clustering))
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
    SG_normal <- single_alt_l1_normalisation(t(currentS[[v]] %*% t(currentG[[v]])))
    #F_normal <- single_alt_l1_normalisation(currentF[[v]])
    #currentF[[v]] <- F_normal$newMatrix
    currentF[[v]] <- currentF[[v]] %*% (SG_normal$Q)
    currentS[[v]] <- solve(SG_normal$Q) %*% currentS[[v]]

    # Update G
    currentG[[v]] <- update_G(Xinput = Xinput[[v]],
                              Finput = currentF[[v]],
                              Sinput = currentS[[v]],
                              Ginput = currentG,
                              psi = psi, k = v)
  
    # Normalise G and S
    FS_normal <- single_alt_l1_normalisation(currentF[[v]] %*% currentS[[v]])
    #G_normal <- single_alt_l1_normalisation(currentG[[v]])
    currentG[[v]] <- currentG[[v]] %*% (FS_normal$Q)
    currentS[[v]] <- currentS[[v]] %*% solve(FS_normal$Q)
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