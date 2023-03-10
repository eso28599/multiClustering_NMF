#generate data
#' In this script, we will create simulations by generating block-diagonal matrices as X matrices
#' 
#' That is instead of generating separately F, S and G and then taking X = FSG^T


# Libraries
library(MASS)
library(Matrix)

## Scratch notes ## 
# Number of row-clusters and column-clusters
n_rowClusters <- 3
n_colClusters <- 3
# Introduce simulated data-views- where we store the views
X_trial <- vector("list", length = 2)

#list to store row clusterings for each dataview
true_row_clusterings <- vector("list", length = length(X_trial)) 
true_col_clusterings <- vector("list", length = length(X_trial))

## Generate first data-view

# Dimensions of each block - how many individuals/features in each cluster
# a vector of length = n_rowClusters
row_dim <- c(100, 200, 300)
col_dim <- c(50, 100, 250)

# Create a list as input for block-diagonal data generation
# list length of no of clusters in this first view
inputList <- vector("list", length = max(length(row_dim), length(col_dim)))
# generate a mvn with n=number of individuals of view and no of features
# equal to number of features of the view
# mean of 5 for each column
# covariance matrix is identity - each feature is independent 
inputList[[1]] <- mvrnorm(n =row_dim[1], mu = rep(5,col_dim[1]), Sigma = diag(col_dim[1]))
inputList[[2]] <- mvrnorm(n =row_dim[2], mu = rep(5,col_dim[2]), Sigma = diag(col_dim[2]))
inputList[[3]] <- mvrnorm(n =row_dim[3], mu = rep(5,col_dim[3]), Sigma = diag(col_dim[3]))

# data for first view - no noise between clusters
X_trial[[1]] <- as.matrix(bdiag(inputList)) 
# generate noise, MVN with no correlation, zero mean and variance of 10
X_noise <- mvrnorm(n = nrow(X_trial[[1]]), mu = rep(0, ncol(X_trial[[1]])), Sigma = diag(10,ncol(X_trial[[1]])))
# add noise to first view
X_trial[[1]] <- X_trial[[1]] + X_noise 

#define vectors indicating true row/column cluster membership for each view
true_row_clusterings[[1]] <- c(rep(1,100),
                               rep(2,200),
                               rep(3,300))
true_col_clusterings[[1]] <- c(rep(1,50),
                               rep(2,100),
                               rep(3,250))

# Repeat for second data-view (with different column dimensions) - same row clustering
# vector of length = n_rowClusters
row_dim <- c(100, 200, 300)
col_dim <- c(100, 100, 100)
# Create a list as input for block-diagonal data generation
inputList <- vector("list", length = max(length(row_dim), length(col_dim)))
inputList[[1]] <- mvrnorm(n =row_dim[1], mu = rep(6,col_dim[1]), Sigma = diag(col_dim[1]))
inputList[[2]] <- mvrnorm(n =row_dim[2], mu = rep(6,col_dim[2]), Sigma = diag(col_dim[2]))
inputList[[3]] <- mvrnorm(n =row_dim[3], mu = rep(6,col_dim[3]), Sigma = diag(col_dim[3]))
X_trial[[2]] <- as.matrix(bdiag(inputList))
X_noise <- mvrnorm(n = nrow(X_trial[[2]]), mu = rep(0, ncol(X_trial[[2]])), Sigma = diag(10,ncol(X_trial[[2]])))
X_trial[[2]] <- X_trial[[2]] + X_noise
true_row_clusterings[[2]] <- c(rep(1,100),
                               rep(2,200),
                               rep(3,300))
true_col_clusterings[[2]] <- c(rep(1,100),
                               rep(2,100),
                               rep(3,100))

# we now have simulated data to run on
#### Run Restrictive Multi-NMTF ####


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

#phi<-Matrix(phi,sparse=TRUE)
#quicker not using sparse matrices!!!

startTime <- Sys.time()
X_trial_NMTF_new <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)
eval_measures_NMTF_new <- evaluate_simulation(X_nmtf = X_trial_NMTF_new,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)
endTime <- Sys.time()
print(endTime-startTime)

#with sparse matrices
#faster with non-sparse 
phi<-Matrix(phi,sparse=TRUE)
psi<-Matrix(psi,sparse=TRUE)
xi<-Matrix(xi,sparse=TRUE)
startTime_sp <- Sys.time()
X_trial_NMTF_new_sparse <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)
eval_measures_NMTF_new_sp <- evaluate_simulation(X_nmtf = X_trial_NMTF_new,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)
endTime_sp <- Sys.time()
print(endTime_sp-startTime_sp)

#try with better normalisation step
#try until convergence
startTime_sp <- Sys.time()
X_trial_NMTF_new_norm <- restMultiNMTF_run(Xinput = X_trial,
                                        Finput = Finit,
                                        Sinput = Sinit,
                                        Ginput = Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi)
eval_measures_NMTF_new_norm <- evaluate_simulation(X_nmtf = X_trial_NMTF_new,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)
endTime_sp <- Sys.time()
print(endTime_sp-startTime_sp)


init_x_trial <- init_mats(X_trial,c(3,3,3),c(3,3,3))
#try with no initialisation or xi
startTimeOld <- Sys.time()
X_trial_NMTF_old <- restMultiNMTF_run(Xinput = X_trial,
                                        KK=c(3,3,3),
                                        LL=c(3,3,3),
                                        phi = phi)
eval_measures_NMTF_old <- evaluate_simulation(X_nmtf = X_trial_NMTF_old,
                                            true_row_clustering = true_row_clusterings,
                                            true_col_clustering = true_col_clusterings)
endTimeOld <- Sys.time()
print(endTimeOld-startTimeOld)


#each dataset with 200 samples, same clusters across views
rowClusters1<-list(c(75,50,30),c(75,50,30),c(75,50,30))
#100, 50,250 features respectively
colClusters1<-list(c(30,30,20),c(10, 20,15),c(15,40,14))


#can change the noise parameter here for level of noise in views
our_data <- multi_view(rowClusters1,colClusters1,noise=5,seed=30)
data_views <- our_data$data_views
init_dv <- init_mats(data_views,c(3,3,3),c(3,3,3))
# Tuning parameters - initialise as all zeros for now
phi <- matrix(0, nrow = length(data_views), ncol = length(data_views))
xi <- matrix(0, nrow = length(data_views), ncol = length(data_views))
psi <- matrix(0, nrow = length(data_views), ncol = length(data_views))
phi[1,2:3]<-0.5

data_views_NMTF <- restMultiNMTF_run(Xinput = data_views,
                                        Finput = init_dv$Finit,
                                        Sinput = init_dv$Sinit,
                                        Ginput = init_dv$Ginit,
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)
eval_dv<- evaluate_simulation(X_nmtf = data_views_NMTF,
                                            true_row_clustering = our_data$truth_rows,
                                            true_col_clustering = our_data$truth_cols)


data_views_NMTF2 <- restMultiNMTF_run(Xinput = data_views,
                                        KK = c(3,3,3),
                                        LL = c(3,3,3),
                                        phi = phi,
                                        xi = xi,
                                        psi = psi,
                                        nIter = nIter)
eval_dv2<- evaluate_simulation(X_nmtf = data_views_NMTF2,
                                            true_row_clustering = our_data$truth_rows,
                                            true_col_clustering = our_data$truth_cols)
