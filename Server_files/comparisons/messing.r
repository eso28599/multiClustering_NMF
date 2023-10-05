library("readxl")
library("rio")
source("sim_studies/increasing_bicl/increasing_bicl/sim_parameters.r")
source("sim_studies/resNMTF_funcs.r")
source("sim_studies/data_generation.r")
import_matrix <- function(filename, col_n = FALSE){
    return(lapply(import_list(filename, col_names=col_n), function(x) as.matrix(x)))
}
#look at comparing posterior distributions
#import as list of matrices instead of as a list
#check it works for overlap and matches
row_cl_dims <- rep(200, 2)
col_cl_dims <- c(100, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 2 biclusters
row_start <- list(c(1, 76), c(1, 76))
row_end <- list(c(75, 200), c(75, 200))
col_start <- list(c(1, 31), c(1, 101))
col_end <- list(c(30, 100), c(100, 250))
row_start_3 <- list(c(1, 76, 151), c(1, 76, 151))
row_end_3 <- list(c(75, 150, 200), c(75, 150, 200))
col_start_3 <- list(c(1, 31, 61), c(1, 101, 151))
col_end_3 <- list(c(30, 60, 100),  c(100, 150, 250))
phi <- matrix(0, nrow = 2, ncol = 2)
phi[1, 2] <- 1
phi[1, c(2,3)] <- 1
noise_level <- 10
i <- 3
check <- one_view_adv(row_dims[i], col_dims[i],
                 row_start[[i]], row_end[[i]],  col_start[[i]],
                  col_end[[i]])

row_cl_dims <- rep(200, 3)

#100, 50,250 features respectively
#col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100), c(40, 50, 60))
col_cl_dims <- c(100, 50, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 2 biclusters
row_start_2 <- list(c(1, 76), c(1, 76), c(1, 76))
row_end_2 <- list(c(75, 200), c(75, 200), c(75, 200))
col_start_2 <- list(c(1, 31), c(1, 11), c(1, 101))
col_end_2 <- list(c(30, 100), c(10, 50), c(100, 250))
# 3 biclusters
row_start_3 <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end_3 <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
col_start_3 <- list(c(1, 31, 61), c(1, 11, 31), c(1, 101, 151))
col_end_3 <- list(c(30, 60, 100), c(10, 30, 50), c(100, 150, 250))
# 4 biclusters
row_start_4 <- list(c(1, 76, 111, 151), c(1, 76, 111, 151), c(1, 76, 111, 151))
row_end_4 <- list(c(75, 110, 150, 200), c(75, 110, 150, 200),
                 c(75, 110, 150, 200))
col_start_4 <- list(c(1, 31, 46, 61), c(1, 11, 21, 31), c(1, 101, 126, 151))
col_end_4 <- list(c(30, 45, 60, 100), c(10, 20, 30, 50), c(100, 125, 150, 250))
image(check$view)
source("sim_studies/returned_data/resNMTF_funcs.r")
data_all <- multi_view(row_cl_dims,
        col_cl_dims, row_start_4, row_end_4, col_start_4, col_end_4, 
        noise=1,row_same_shuffle=TRUE, col_same_shuffle=rep(FALSE,3), seed=5)
data_views <- data_all$data_views
phi_mat <- matrix(0, 3, 3)
phi_mat[1,c(1,2)] <- 10000
phi_mat[2,c(3)] <- 10000
t1 <- Sys.time()
nmtf_results <- restMultiNMTF_run(Xinput = data_views,
         phi = phi_mat, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))
check_score <- new_sil_score(data_views[1:2], nmtf_results$Foutput,
     nmtf_results$Goutput, nmtf_results$row_clusters, nmtf_results$col_clusters, 2)
jaccard_results(nmtf_results$row_clusters[[1]],nmtf_results$col_clusters[[1]],
         data_all$truth_rows[[1]], data_all$truth_cols[[1]])
nmtf_results$Soutput

mage(data_views[[3]])
#return plots of images
i <- 1
image(data_views[[3]][nmtf_results$row_clusters[[3]][,i]==1,
        nmtf_results$col_clusters[[3]][,i]==1])
data_views <- import_matrix("sim_studies/returned_data/data.xlsx", col_n=TRUE)
#row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
row_cl_dims <- rep(200, 2)
col_cl_dims <- c(100, 250)

row_same_shuffle <- TRUE
col_same_shuffle <- FALSE

# 2 biclusters
row_start_2 <- list(c(1, 76), c(1, 76))
row_end_2 <- list(c(1, 76), c(1, 76))
col_start_2 <- list(c(1, 31), c(1, 101))
col_end_2 <- list(c(30, 100), c(100, 250))


check <- one_view_adv(200, 100,
                 c(1, 76), c(75, 200),
                   c(1, 31),
                  c(30, 100))
image(check[[1]])
image(as.matrix(dist(check[[1]])))
row_dims <-  rep(200, 3)
#100, 100,250 features respectively
col_dims <- c(100, 100, 250)

row_start <- list(c(1, 76, 151), c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 200), c(75, 150, 200), c(75, 150, 200))
col_start <- list(c(1, 31, 61), c(1, 31, 61), c(1, 101, 151))
col_end <- list(c(30, 60, 100),c(30, 60, 100), c(100, 150, 250))


check <- one_view_adv(200, 100,
                 c(1, 76, 151), c(75, 150, 200),
                   c(1, 31, 61),
                  c(30, 60, 100))


source("sim_studies/increasing_bicl/resNMTF_funcs.r")
phi_mat <- matrix(0, 2, 2)
phi_mat[1,2] <- 500000
t1 <- Sys.time()
nmtf_results <- restMultiNMTF_run(Xinput = data_views[1:2], 
         phi = phi_mat, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))

check_score <- new_sil_score(data_views[1:2], nmtf_results$Foutput,
     nmtf_results$Goutput, nmtf_results$row_clusters, nmtf_results$col_clusters, 2)


#generate data 



n_views <- length(data_views)
n_clusts <- dim(nmtf_results$Foutput[[1]])[2]
data_messed <- vector(mode = "list", length = n_views)
for(i in 1:n_views){
      #correct shuffling
    dims <- dim(data_views[[i]])
    data_messed[[i]] <- matrix(sample(data_views[[i]]),
                              dims[1], dims[2])
    data_messed[[i]] <- apply(data_messed[[i]], 2, function(x) x/sum(x))
}

t1 <- Sys.time()
nmtf_mess <- restMultiNMTF_run(Xinput = data_messed, 
         phi = phi_mat, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))




row_cl_dims <- list(c(75, 75, 50, 50), c(75, 75, 50, 50), c(75, 75, 50, 50))
#100, 50,250 features respectively
col_cl_dims <- list(c(30, 30, 40, 20), c(10, 20, 20, 30), c(100, 50, 100, 30))

row_cl_dims <- list(c(75, 75, 50), c(75, 75, 50), c(75, 75, 50))
#100, 50,250 features respectively
col_cl_dims <- list(c(30, 30, 40), c(10, 20, 20), c(100, 50, 100))

row_cl_dims <- list(c(75, 50, 50, 25), c(75, 50, 50, 25))
#100, 50,250 features respectively
col_cl_dims <- list(c(10, 20, 10, 10), c(100, 50, 75, 25))
noise_level <- 3
n_views <- length(row_cl_dims)
phi <- matrix(0, nrow = n_views, ncol = n_views)
phi[1, 2] <- 1
phi_mat <- 0.01*phi
#read in data
data_all <- multi_view(row_cl_dims,
        col_cl_dims,noise=noise_level,row_same_shuffle=TRUE, col_same_shuffle=FALSE, seed=10)
data_views <- data_all$data_views
data_sorted <- vector(mode = "list", length = length(data_views))
for(i in 1:length(data_views)){
    #correct shuffling
    data_sorted[[i]] <- apply(data_views[[i]], 2, function(x) x/sum(x))

}
#synthetic data results
t1 <- Sys.time()
nmtf_results <- restMultiNMTF_run(Xinput = data_views, 
         phi = phi_mat, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))

t1 <- Sys.time()
nmtf_results_mess <- restMultiNMTF_run(Xinput = data_messed, 
         phi = phi_mat, relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))

evaluate_simulation_comp(nmtf_results$row_clusters, nmtf_results$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_sorted, index = 2)

jaccard_results(nmtf_results$row_clusters[[1]],nmtf_results$col_clusters[[1]],
         data_all$truth_rows[[1]], data_all$truth_cols[[1]])

#new clustering results - soft clustering 
nmtf_results$Sil_score

t1 <- Sys.time()
test_check <- check_biclusters(data_sorted, nmtf_results$Foutput, 5, 2)
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))

t1 <- Sys.time()
clust_res <- clustering_res_NMTF2(data_sorted, nmtf_results$Foutput,
                        nmtf_results$Goutput, nmtf_results$Soutput, 5, 2)
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))



t1 <- Sys.time()
find_thres <- get_thresholds(data_sorted, nmtf_results$Foutput, 5, 2)
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))

find_thres[[1]][, 1]
test_scores <- new_sil_score(data_sorted,nmtf_results$Foutput,
     clust_res$row_clustering, clust_res$col_clustering, 2)

#check for diff number of clusters
k_input <- 4 * rep(1, length = n_views)
testing <- restMultiNMTF_run(Xinput = data_sorted, phi = phi_mat,
                    KK = k_input,  LL = k_input, relErr = "Sil")
clust_four <- clustering_res_NMTF2(testing$Foutput, 
                        testing$Goutput, testing$Soutput, 2)
four_scores <- new_sil_score(data_sorted,testing$Foutput,
     clust_four$row_clustering, clust_four$col_clustering, 2)

indices <- clust_res$row_clustering[[1]][, 1] == 1
not_ind <- clust_res$row_clustering[[1]][, 1] != 1
try <- dist(data_sorted[[1]][c((1:200)[indices], (1:200)[not_ind]),
            clust_res$col_clustering[[1]][,1] == 1], diag = T, upper = F)
heatmap(as.matrix(try))
idx <- (1:200)[data_all$truth_rows[[1]]==3]
idx_n <- (1:200)[data_all$truth_rows[[1]]!=3]
heatmap(data_sorted[[1]][c(idx),
            clust_res$col_clustering[[1]][,1] == 1])

as.matrix(as.dist(cor(t(data_sorted[[1]][1:20,1:10]), method = "spearman"), upper = TRUE, diag = TRUE))[1:5,1:5]
n_views <- length(data_sorted)
n_clusts <- dim(nmtf_results$Foutput[[1]])[2]
data_messed <- vector(mode = "list", length = n_views)
data_mvn <- vector(mode = "list", length = n_views)
for(i in 1:n_views){
      #correct shuffling
    dims <- dim(data_sorted[[i]])
    data_messed[[i]] <- matrix(sample(data_sorted[[i]]),
                              dims[1], dims[2])
    data_messed[[i]] <- apply(data_messed[[i]], 2, function(x) x/sum(x))
    data_mvn[[i]] <- abs(mvrnorm(dims[1],
           mu = colMeans(data_sorted[[i]]), Sigma = cov(data_sorted[[i]])))
    data_mvn[[i]] <- apply(data_mvn[[i]], 2, function(x) x/sum(x))
}

n_views <- length(data_sorted)
n_clusts <- dim(nmtf_results$Foutput[[1]])[2]
k_input <- n_clusts * rep(1, length = n_views)
testing <- restMultiNMTF_run(Xinput = data_messed, phi = phi_mat,
                    KK = k_input,  LL = k_input, relErr = "Sil")
testing_mvn <- restMultiNMTF_run(Xinput = data_mvn, phi = phi_mat,
            KK = k_input,  LL = k_input, relErr = "Sil")
clust_mvn <- clustering_res_NMTF2(testing_mvn$Foutput, 
                        testing_mvn$Goutput, testing_mvn$Soutput, 2)
x <- nmtf_results$Foutput[[1]][, 1]
x1_noise <- testing$Foutput[[1]][, 1]
x2_noise <- testing_mvn$Foutput[[1]][, 1]

KL(rbind(x1_noise, x), unit = "log2")
JSD(rbind(x,x1_noise), unit = "log2")
JSD(rbind(x2_noise,x1_noise), unit = "log2")

testing2 <- restMultiNMTF_run(Xinput = data_messed, phi = phi_mat, relErr = "Sil")
testing_mvn2 <- restMultiNMTF_run(Xinput = data_mvn, phi = phi_mat, relErr = "Sil")


x1_noise2 <- testing2$Foutput[[1]][, 1]
x2_noise2 <- testing_mvn2$Foutput[[1]][, 1]
KL(rbind(x1_noise2, x), unit = "log2")
JSD(rbind(x,x1_noise2), unit = "log2")
JSD(rbind(x,x2_noise2), unit = "log2")
JSD(rbind(x1_noise2,x1_noise), unit = "log2")
JSD(rbind(x2_noise,x1_noise2), unit = "log2")
JSD(rbind(x2_noise,x2_noise2), unit = "log2")

x1_noise3 <- testing$Foutput[[1]][, 1]
x2_noise2 <- testing_mvn2$Foutput[[1]][, 1]
JSD(rbind(x1_noise3, x), unit = "log2")
JSD(rbind(x,x1_noise2), unit = "log2")
JSD(rbind(x,x2_noise2), unit = "log2")
JSD(rbind(x1_noise2,x1_noise), unit = "log2")
JSD(rbind(x1_noise3,x1_noise), unit = "log2")
JSD(rbind(x1_noise2,x1_noise3), unit = "log2")
JSD(rbind(x2_noise,x1_noise2), unit = "log2")
JSD(rbind(x2_noise,x2_noise2), unit = "log2")


data_views <- import_matrix("comparisons/scen_1/data/data.xlsx", col_n=TRUE)
true_rows <- import_matrix("comparisons/scen_1/data/data.xlsx", col_n=TRUE)
true_cols <- import_matrix("comparisons/scen_1/data/data.xlsx", col_n=TRUE)
#present results
data_sorted <- vector(mode = "list", length = length(data_views))
for(i in 1:length(data_views)){
    #correct shuffling
    data_sorted[[i]] <- apply(data_views[[i]], 2, function(x) x/sum(x))

}

nmtf_results2 <- restMultiNMTF_run(Xinput = data_sorted, phi = phi_mat, relErr = "Sil")
clust_res <- clustering_res_NMTF(nmtf_results2$Foutput, 
                        nmtf_results2$Goutput, nmtf_results2$Soutput, 2)

data_messed<- vector(mode = "list", length = length(data_views))
for(i in 1:length(data_views)){
    #correct shuffling
    dims <- dim(data_views[[i]])
    data_messed[[i]] <- matrix(sample(apply(data_views[[i]], 2, function(x) x/sum(x))),
                            dims[1], dims[2])

}

nmtf_results2$
colSums(data_sorted[[1]])
heatmap(data_sorted[[1]])
a <- matrix(c(4,5,3,4,2,8,1,2,3), nrow=3)
sort_rels <- function(mat){
   sorted <- rep(0, n=3)
    for(i in 1:3){
        sort_links <- rep(0, n=3)
        for(j in 1:3){
            sort_links[j] <- sum(abs(mat[i,] - mat[i,j]))
        }
        sorted[i] <- which.max(sort_links)
    } 
    return(sorted)
}
sort_rels(nmtf_results$Soutput[[2]])
a1 <- sum(abs(a[1,] - a[1,1]))
colSums( t(nmtf_results2$Goutput[[1]])) - rowSums(data_sorted[[1]])

SG_normal <- single_alt_l1_normalisation(t(nmtf_results2$Soutput[[1]] %*% t(nmtf_results2$Goutput[[1]])))
FS_normal <- single_alt_l1_normalisation(nmtf_results2$Foutput[[1]] %*% nmtf_results2$Soutput[[1]])

nmtf_results <- restMultiNMTF_run(Xinput = data_views, phi = phi_mat, relErr = "Sil")
qheatmap(nmtf_results$Foutput[[1]]%*%nmtf_results$Soutput[[1]]%*%t(nmtf_results$Goutput[[1]]))
heatmap(nmtf_results2$Foutput[[1]]%*%nmtf_results2$Soutput[[1]]%*%t(nmtf_results2$Goutput[[1]]))

nmtf_results_messed <- restMultiNMTF_run(Xinput = data_messed, phi = phi_mat, relErr = "Sil")


plot((nmtf_results_messed$Foutput[[1]]/colSums(nmtf_results_messed$Foutput[[1]]))[,1])
rowSums(nmtf_results2$Foutput[[1]]%*%nmtf_results2$Soutput[[1]]%*%diag(colSums(nmtf_results2$Goutput[[1]])))
max(rowSums(nmtf_results$Foutput[[1]]%*%nmtf_results$Soutput[[3]]))
colSums(nmtf_results$Foutput[[1]]%*%nmtf_results$Soutput[[1]])
x <- abs((nmtf_results_messed$Foutput[[1]]/colSums(abs(nmtf_results_messed$Foutput[[1]])))[,1])
x_2 <- (nmtf_results_messed$Foutput[[1]]/colSums(nmtf_results_messed$Foutput[[1]]))[,2]
x2 <- (nmtf_results$Foutput[[1]]/colSums(nmtf_results$Foutput[[1]]))[,1]
x3 <- (nmtf_results$Foutput[[1]]/colSums(nmtf_results$Foutput[[1]]))[,2]

KL_score <- rep(0,100)
for(i in 1:100){
    x_ex <- rnorm(200, mean(x) ,sd(x))
    x_ex <- abs(x_ex)/sum(abs(x_ex))
    KL_score[i] <- KL(rbind(x_ex, x/sum(x)))
}

plot( (abs(x3)/(sum(abs(x3))))[abs(x3)/(sum(abs(x3)))<0.005])

sum(apply(nmtf_results2$Soutput[[1]]%*%t(nmtf_results2$Goutput[[1]]), 2, which.max)==3)
rowSums(nmtf_results$Soutput[[1]]%*%t(nmtf_results$Goutput[[1]]))
colSums(nmtf_results2$Soutput[[3]]%*%t(nmtf_results2$Goutput[[3]]))
apply(nmtf_results2$Soutput[[1]]%*%t(nmtf_results2$Goutput[[1]]), 2, which.max)
sum(nmtf_results2$Foutput[[1]]> 0.001)
qtest1 <- nmtf_results$Foutput[[2]] %*% nmtf_results$Soutput[[2]] %*% t(nmtf_results$Goutput[[2]])
n <- 3
G <- apply(nmtf_results$Goutput[[n]]/apply(nmtf_results$Goutput[[n]],1,max) ==1, 2, as.numeric)
F <- apply(nmtf_results$Foutput[[n]]/apply(nmtf_results$Foutput[[n]],1,max) ==1, 2, as.numeric)
S <- nmtf_results$Soutput[[n]]
test_mat <- F %*% S %*% t(G)
colSums(nmtf_results2$Soutput[[3]]%*%t(nmtf_results2$Goutput[[3]]))
test_mat <- clustering_res_NMTF_update(nmtf_results$Foutput,
             nmtf_results$Goutput, 
             nmtf_results$Soutput, n_views, 0.5)

openxlsx::write.xlsx(test_mat$row_clustering, file = "path_to_save.xlsx")

test_mat1 <- clustering_res_NMTF(nmtf_results$Foutput,
             nmtf_results$Goutput, 
             nmtf_results$Soutput, n_views)

install.packages("/Users/ellaorme/Downloads/mvcluster_1.0.tar.gz", repos=NULL, type='source')
install.packages("/mvcluster_1.0.tar.gz", repos=NULL, type='source')
library(devtools)
devtools::install_local("/Users/ellaorme/Downloads/mvcluster_1.0.tar.gz")

devtools::install_local("/Users/ellaorme/Downloads/RcppArmadillo_0.11.2.4.0.tar.gz")
install.packages("/Users/ellaorme/Downloads/RcppArmadillo_0.11.2.4.0.tar.gz", repos=NULL, type='source')

install_version("RcppArmadillo", version = "0.11.2.4.0", repos = "http://cran.us.r-project.org")
library(GFA)
norm <- normalizeData(data_views,type="scaleFeatures")
## Get the model options - we want to detect bicluster structure
opts <- getDefaultOpts(bicluster=TRUE)
opts$convergenceCheck <- TRUE #We want to check for sampling chain convergence
## Giving a vague prior stating that roughly 1/4 of the variance is noise
opts <- informativeNoisePrior(data_views,opts,noiseProportion=0.25,conf=0.1)

## Infer the model

res2 <- gfa(norm$train,opts=opts,K=8) # supress messages
res2$W #seems to be coefficients membership of col clusters? although negative don;t sum to one etc
res2$Z #denotes col cluster membership but all variables appear to belong to W



row_dims <- 80
col_dims <- 100
row_start <- c(1, 41, 55)
row_end <- c(40, 60, 70)
col_start <- c(1, 51, 75)
col_end <- c(40, 80, 90)

check_mat <- one_view_adv(row_dims, col_dims, row_start, row_end,  col_start, col_end, noise =2)
image(check_mat$view)
original <- one_view(c(40,30,15), c(40,40,10),noise=10)

#test one view
#synthetic data result
phi_mat <- matrix(0, 1, 1)
t1 <- Sys.time()
nmtf_results_ov <- restMultiNMTF_run(Xinput = list(check_mat$view), 
         phi = matrix(0,1,1), relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))

t1 <- Sys.time()
nmtf_results_og <- restMultiNMTF_run(Xinput = list(original$view), 
         phi = matrix(0,1,1), relErr = "Sil")
t2 <- Sys.time()
print(paste0("time taken was ", t2-t1))


#check it works for overlap and matches
row_dims <- 80
col_dims <- 100
row_start <- c(1, 41, 55)
row_end <- c(40, 60, 70)
col_start <- c(1, 51, 75)
col_end <- c(40, 80, 90)