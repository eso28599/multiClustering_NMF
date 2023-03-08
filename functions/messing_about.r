#getting running and figuring out why update steps affect

attempt_old<- update_F(X_trial[[1]], Finit, Sinit[[1]], Ginit[[1]], phi, 1)
attempt_new<- update_F1(X_trial[[1]], Finit, Sinit[[1]], Ginit[[1]], phi, 1)

#checking if difference is introduced by normalisation - it isn't 
Q_F_old <- single_alt_l1_normalisation(attempt_old)$Q
Q_F_new <- single_alt_l1_normalisation(attempt_new)$Q
normal_old <- attempt_old %*% solve(Q_F_old)
normal_new <- attempt_new %*% solve(Q_F_new)
paste("total absoulte difference between updates is", sum(abs(normal_old-normal_new)))
paste("max relative difference between updates is",max(abs(normal_old-normal_new)))
#check with different values
attempt_old<- update_F(X_trial[[1]], X_trial_NMTF$Foutput, X_trial_NMTF$Soutput[[1]],X_trial_NMTF$Goutput[[1]], phi, 1)
attempt_new<- update_F1(X_trial[[1]], X_trial_NMTF$Foutput, X_trial_NMTF$Soutput[[1]],X_trial_NMTF$Goutput[[1]], phi, 1)

paste("total absoulte difference between updates is", sum(abs(attempt_old-attempt_new)))
paste("max relative difference between updates is",max(abs(attempt_old-attempt_new)/abs(attempt_old)))

#GS^T - calculate once to improve speed
com_mat1 <- Ginit[[1]] %*% t(Sinit[[1]])
numerator_matrix1 <- X_trial[[1]] %*% com_mat1
denominator_matrix1 <- Finit[[1]] %*% t(com_mat1) %*% com_mat1

#calculate the column vector based on phi that is needed
phi_vec1 <- (phi+t(phi))[,1]
num_mat_prod1 <- star_prod(phi_vec1, Finit)
denom_mat_prod1 <- sum(phi_vec1)* Finit[[1]]

outputF <- Finit[[1]]* ((numerator_matrix1 + num_mat_prod1) / (denominator_matrix1 + denom_mat_prod1))

com_mat1 <- X_trial_NMTF$Goutput[[1]] %*% t(X_trial_NMTF$Soutput[[1]])
numerator_matrix1 <- X_trial[[1]] %*% com_mat1
denominator_matrix1 <- X_trial_NMTF$Foutput[[1]] %*% t(com_mat1) %*% com_mat1

#calculate the column vector based on phi that is needed
phi_vec1 <- (phi+t(phi))[,1]
num_mat_prod1 <- star_prod(phi_vec1, X_trial_NMTF$Foutput)
denom_mat_prod1 <- sum(phi_vec1)* X_trial_NMTF$Foutput[[1]]
denominator_matrix1 + denom_mat_prod1

outputF <- Finit[[1]]* ((numerator_matrix1 + num_mat_prod1) / (denominator_matrix1 + denom_mat_prod1))




#try and get it running
X_trial_NMTF <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)


startTime <- Sys.time()
X_trial_NMTF$
#check the difference in errors

X_trial_NMTF_new <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)

endTime <- Sys.time()

print(endTime-startTime)

startTime2 <- Sys.time()

X_trial_NMTF_old <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)

endTime2 <- Sys.time()
print(endTime2-startTime2)

paste("total absoulte difference between updates is", sum(abs(X_trial_NMTF_old$Foutput[[1]]-X_trial_NMTF_new$Foutput[[1]])))
paste("max difference between updates is",max(abs(X_trial_NMTF_old$Error-X_trial_NMTF_new$Error)))
#also with G changed
startTime <- Sys.time()

X_trial_NMTF_new <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)

endTime <- Sys.time()

X_trial_NMTF_new$Foutput
apply(X_trial_NMTF_new$Foutput[[1]], 1, which.max)

eval_measures_NMTF_new <- evaluate_simulation(X_nmtf = X_trial_NMTF_new,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)

eval_measures_NMTF_old <- evaluate_simulation(X_nmtf = X_trial_NMTF_old,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)

#check if re-labelling affects - double check but doesn't seem to
# as long as solve_LSAP function always works properly
true_row_clusterings2<-ifelse(true_row_clusterings[[1]]==1,"a",ifelse(true_row_clusterings==2,"b","c"))
true_row_clusterings2<-ifelse(true_row_clusterings2=="a",3,ifelse(true_row_clusterings2=="b"),2,1)
eval_measures_NMTF_old2 <- evaluate_simulation(X_nmtf = X_trial_NMTF_old,
                                            true_row_clustering = true_row_clusterings2,
                                            true_col_clustering = true_col_clusterings)


#they are going to different local minima
#smaller errors for the solution achieved by theo's method
#is this just for this one run?
#need to do a simulation study for this
#is it just accumulation of small numerical error? This would suggest the method really isn't robust at all 
#see if with same data but two different sets of small noise we get the same results - with either method

# if not where is it being introduced?
# do we get the same results if we use the simulated data based off noisy FSG^T construction?

# run the different functions at different starting points - is it robust

X_trial <- data_views
# Initialisation of F, S, G lists - one matrix for each data view
Finit <- vector("list", length = length(X_trial))
Sinit <- vector("list", length = length(X_trial))
Ginit <- vector("list", length = length(X_trial))

#initialise based on svd - why?
for (i in 1:length(X_trial)){
  ss <- svd(X_trial[[i]])
  Finit[[i]] <- ss$u
  Sinit[[i]] <- diag(ss$d)
  Ginit[[i]] <- ss$v
}

#numbers of clusters in each view - doesn't need the last element I don't think 
KK <- c(3,3,3)
LL <- c(3,3,3)
#changing initialisation of S to be radnom 
for (i in 1:length(X_trial)){
  L <- LL[[i]]
  K <- KK[[i]]
  Sinit[[i]] <- mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L))
}
# Sinit is no longer strictly positive? 
# Select top K for row-clusters and L for column clusters

for (i in 1:length(X_trial)){
  L <- LL[[i]]
  K <- KK[[i]]
  Finit[[i]] <- abs(Finit[[i]][,1:K])
  Sinit[[i]] <- abs(Sinit[[i]][1:K,1:L])
  Ginit[[i]] <- abs(Ginit[[i]][,1:L])
}

# Tuning parameters - initialise as all zeros for now
phi <- matrix(0, nrow = length(X_trial), ncol = length(X_trial))
xi <- matrix(0, nrow = length(X_trial), ncol = length(X_trial))
psi <- matrix(0, nrow = length(X_trial), ncol = length(X_trial))
nIter <- 1000
phi[1,2] <- 0.5
phi[1,3] <- 0.5
phi<-Matrix(phi,sparse=TRUE)
psi<-Matrix(psi,sparse=TRUE)


new_parameters1 <- restrictiveMultiNMTF_algo(X = data_views, 
                                                  R = 1, 
                                                  Finit = Finit, 
                                                  Sinit = Ginit,
                                                  Ginit = Sinit,
                                                  phi = phi,
                                                  xi = xi, 
                                                  psi = psi)  
data_views_NMTF_newFGS <- restrictiveMultiNMTF(Xinput = data_views,
                                        R=1,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)


data_views_NMTF_newFGS <- restrictiveMultiNMTF(data_views,
                                        R=1,
                                         Finit,
                                       Sinit,
                                      Ginit,
                                         phi,
                                         xi,
                                        psi,
                                     nIter)                                   
                    

star_prod(psi_vec,Ginit)
