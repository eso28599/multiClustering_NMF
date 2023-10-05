#' An example of running Multi-view restrictive Non-negative Matrix Tri-Factorisation
#' 
#' Use dataExample as input matrix
#' 

# Directory
setwd("/Users/ellaorme/GitHub/Data_integration")
source("/Users/ellaorme/GitHub/Data_integration/multiClustering_NMF/functions/restrictive_NMTF_singularCheck.R")
# Libraries
library(data.table)
library(MASS)
library(clue)
library(fossil)
library(mclust)
library(aricode)
library(whitening)

# Data
# dataExample: List 
# labelsExample
#read in data set 2b 
# only consider expression and mirna data 
#dataExample <- combined_set_2b[c(1,3)]
##labelsExample <- combined_set_2b[4][[1]]
dataExample <- data
labelsExample <- subtypes
#### RUN Restrictive Multi-NMTF ####
# True clustering

true_row_clusterings <- vector("list", length = length(dataExample))
#true_col_clusterings <- vector("list", length = length(dataExample))
names(true_row_clusterings) <- names(dataExample)
for (i in 1:length(dataExample)){
  true_row_clusterings[[i]] <- labelsExample
  #true_col_clusterings[[i]] <- rep(1, ncol(dataExample[[i]]))
}

# Initialization
Finit <- vector("list", length = length(dataExample))
Sinit <- vector("list", length = length(dataExample))
Ginit <- vector("list", length = length(dataExample))
#initialise matricies by svd 
for (i in 1:length(dataExample)){
  ss <- svd(dataExample[[i]])
  Finit[[i]] <- ss$u
  Sinit[[i]] <- diag(ss$d)
  Ginit[[i]] <- ss$v
}
KK <- rep(length(unique(labelsExample)),length(dataExample)) #vector of no of clusters in each of the omics
LL <- rep(2, length(dataExample)) #? why 2 #number of column clusters in each omics data,
#is this 2 in the example used? How to change this for ours when this is unknown 
#loop over each OMIC
for (i in 1:length(dataExample)){
  L <- LL[[i]] #this is the same for each i?
  K <- KK[[i]]
  Sinit[[i]] <- mvrnorm(n = K, mu = rep(0,L), Sigma = diag(L))
}
# Select top K for row-clusters and L for column clusters


for (i in 1:length(dataExample)){
  L <- LL[[i]]
  K <- KK[[i]]
  Finit[[i]] <- abs(Finit[[i]][,1:K]) #select first K columns
  Sinit[[i]] <- abs(Sinit[[i]][1:K,1:L])
  Ginit[[i]] <- abs(Ginit[[i]][,1:L]) #select first L columns
}

# Tuning parameters # the matrices implementing the restriction matrices
#these are the R matrices indicating which restrictions are in place 
phi <- matrix(0, nrow = length(dataExample), ncol = length(dataExample))
xi <- matrix(0, nrow = length(dataExample), ncol = length(dataExample))
psi <- matrix(0, nrow = length(dataExample), ncol = length(dataExample))
nIter <- 1000

test_phi <- c(0, 0.01, 0.05,0.1, 0.2, 0.5, 0.65, 0.8, 1.0, 5.0, 10.0, 20.0, 50.0)
test_phi <- c(0.4)
tuningAccuracy <- vector("list", length = length(test_phi))
tuningAdjRandValue <- vector("list", length = length(test_phi))
tuningNMIvalue <- vector("list", length = length(test_phi))
tuningError <- vector("list", length = length(test_phi))

#implementation of different parameter values 
#-----------------------------------------
# Assume n_v = 2 data-views
k <- 1
# this imposes the 1st and 2nd data clusters having the same row clusters
#leave as is for now, then will try and change to enforce same row clustering for all views
#we use emp
for (k in 1:length(test_phi)){
  phi[1,2] <- test_phi[k]
  dataExample_NMTF <- restMultiNMTF_run(Xinput = dataExample,
                                            Finput = Finit,
                                            Sinput = Sinit,
                                            Ginput = Ginit,
                                            phi = phi,
                                            xi = xi,
                                            psi = psi,
                                            nIter = nIter)
  eval_measures_dataExample_NMTF <- evaluate_simulation(X_nmtf = dataExample_NMTF,
                                                 true_row_clustering = true_row_clusterings,
                                                 true_col_clustering = true_col_clusterings)
  tuningAccuracy[[k]] <- eval_measures_dataExample_NMTF$accuracy
  tuningAdjRandValue[[k]] <- eval_measures_dataExample_NMTF$ARI
  tuningNMIvalue[[k]] <- eval_measures_dataExample_NMTF$NMI
  tuningError[[k]] <- dataExample_NMTF$Error
  if (k == 1){
    finalAccuracy_single <- tuningAccuracy[[k]]
    finalAdjRandValue_single <- tuningAdjRandValue[[k]]
    finalNMIvalue_single <- tuningNMIvalue[[k]]
    finalError_single <- tuningError[[k]]
  }
  if (test_phi[k] == 1){
    finalAccuracy_allOnes <- tuningAccuracy[[k]]
    finalAdjRandValue_allOnes <- tuningAdjRandValue[[k]]
    finalNMIvalue_allOnes <- tuningNMIvalue[[k]]
    finalError_allOnes <- tuningError[[k]]
  }
}
# Find optimal parameter:
meanACC <- lapply(tuningAccuracy, rowMeans)
chosen_metric <- 1
temp_mean <- meanACC[[1]][1]
for (i in 2:length(tuningAccuracy)){
  if (meanACC[[i]][1] > temp_mean){
    temp_mean <- meanACC[[i]][1]
    chosen_metric <- i
  }
}

phi[1,2] <- test_phi[chosen_metric]
# ----------------------------------------------------------
# basic implementation
phi[1,2] <- 0.5
dataExample_NMTF <- restMultiNMTF_run_noScale(Xinput = dataExample,
                                          Finput = Finit,
                                          Sinput = Sinit,
                                          Ginput = Ginit,
                                          phi = phi,
                                          xi = xi,
                                          psi = psi,
                                          nIter = nIter)


#not working 
test_X <- dataExample[[2]]
check <- colSums(test_X)
problem_point <-  single_alt_l1_normalisation(test_X)

eval_measures_dataExample_NMTF <- evaluate_simulation(X_nmtf = dataExample_NMTF,
                                                      true_row_clustering = true_row_clusterings,
                                                      true_col_clustering = true_col_clusterings)

finalAccuracy_multiple <- eval_measures_dataExample_NMTF$accuracy
finalAdjRandValue_multiple <- eval_measures_dataExample_NMTF$ARI
finalNMIvalue_multiple <- eval_measures_dataExample_NMTF$NMI
finalError_multiple <- dataExample_NMTF$Error

## HSIC
hsic_values <- hsic_independence_test(dataExample)
phi <- hsic_values
phi[lower.tri(phi)] <- 0
diag(phi) <- 0
dataExample_NMTF <- restMultiNMTF_run(Xinput = dataExample,
                                          Finput = Finit,
                                          Sinput = Sinit,
                                          Ginput = Ginit,
                                          phi = phi,
                                          xi = xi,
                                          psi = psi,
                                          nIter = nIter)
eval_measures_dataExample_NMTF <- evaluate_simulation(X_nmtf = dataExample_NMTF,
                                                      true_row_clustering = true_row_clusterings,
                                                      true_col_clustering = true_col_clusterings)
finalAccuracy_hsic <- eval_measures_dataExample_NMTF$accuracy
finalAdjRandValue_hsic <- eval_measures_dataExample_NMTF$ARI
finalNMIvalue_hsic <- eval_measures_dataExample_NMTF$NMI
finalError_hsic <- dataExample_NMTF$Error
