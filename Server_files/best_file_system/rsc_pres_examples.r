#examples for presentation
source("best_file_system/data_generation_toy.r")
source("best_file_system/evaluation_funcs.r")
source("best_file_system/resNMTF_funcs.r")
#toy examples
#a -  same biclusters in terms of columns across two of the views
##view 2 column dimensions becomes 100: (30, 30, 40)#simulation parameters
noise_level <- 5
#each dataset with 200 samples, same clusters across views
row_dims <-  rep(200, 2)
#100, 100,250 features respectively
col_dims <- c(100, 250)

row_start <- list(c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 200),  c(75, 150, 200))
col_start <- list(c(1, 31, 61),  c(1, 101, 151))
col_end <- list(c(30, 60, 100), c(100, 150, 250))
niter <- 1000
# exhausitive, exclusive clusters
n_views <- length(row_dims)
no_row_cl <- sapply(row_dims, length)
no_col_cl <- sapply(col_dims, length)
#parameters for restrictive NMTF
phi <- matrix(0, nrow = n_views, ncol = n_views)
res_val <- 5000
phi[1, 2] <- 1
phi_mat <- res_val * phi

#generate data - return images
data_all <- multi_view(row_dims,
        col_dims, row_start, row_end, col_start, col_end, 
        noise=noise_level,row_same_shuffle=TRUE, col_same_shuffle=c(F, F), seed=10)
source("best_file_system/resNMTF_funcs.r")
nmtf_results <- restMultiNMTF_run(Xinput = data_all$data_views, 
            phi = phi_mat, relErr = "Sil")
#evaluate results
source("best_file_system/evaluation_funcs.r")
nmtf_scores <- evaluate_simulation_comp(nmtf_results$row_clusters , nmtf_results$col_clusters,
 data_all$truth_rows, data_all$truth_cols, data_all$data_views, index = 2)

#produce first image
#reorder rows
for(j in 1:n_views){   
    row_order <-c()
    col_order <- c()

    for(i in 1:length(row_start[[1]])){
        row_order <- c(row_order, 
                (1:((row_end[[j]][3])))[nmtf_results$row_clusters[[j]][, i]==1])
        col_order <- c(col_order, 
                (1:((col_end[[j]][3])))[nmtf_results$col_clusters[[j]][, i]==1])
    }
    image(data_all$data_views[[j]][row_order,col_order])
}
image(data_all$data_views[[1]])
image(data_all$data_views[[2]])

#now have overlapping biclusters
row_start <- list(c(1, 76, 151), c(1, 76, 151))
row_end <- list(c(75, 150, 180),  c(75, 150, 180))
col_start <- list(c(1, 31, 61),  c(1, 101, 151))
col_end <- list(c(40, 60, 85), c(120, 150, 200))
data_overlap <- multi_view(row_dims,
        col_dims, row_start, row_end, col_start, col_end, 
        noise=noise_level,row_same_shuffle=TRUE, col_same_shuffle=c(F, F), seed=10)
source("best_file_system/resNMTF_funcs.r")
nmtf_overlap <- restMultiNMTF_run(Xinput = data_overlap$data_views, 
            phi = phi_mat, relErr = "Sil")
#evaluate results
source("best_file_system/evaluation_funcs.r")
nmtf_scores <- evaluate_simulation_comp(nmtf_overlap$row_clusters , nmtf_overlap$col_clusters,
 data_overlap$truth_rows, data_overlap$truth_cols, data_overlap$data_views, index = 2)

#produce first image
#reorder rows
for(j in 1:n_views){   
    row_order <-c()
    col_order <- c()

    for(i in 1:length(row_start[[1]])){
        row_order <- c(row_order, 
                (1:((row_end[[j]][3])))[nmtf_overlap$row_clusters[[j]][, i]==1])
        col_order <- c(col_order, 
                (1:((col_end[[j]][3])))[nmtf_overlap$col_clusters[[j]][, i]==1])
    }
    image(data_overlap$data_views[[j]][row_order,col_order])
    image(data_overlap$data_views[[j]])
}
image(data_overlap$data_views[[j]])

#view 1
#rows
a <- (1:200)[(nmtf_overlap$row_clusters[[1]][, 1])==1]
b <- (1:200)[(nmtf_overlap$row_clusters[[1]][, 2])==1]
c <- (1:200)[(nmtf_overlap$row_clusters[[1]][, 3])==1]
d <- (1:200)[(nmtf_overlap$row_clusters[[1]][, 4])==1]
no_row <-(1:200)[-sort(union(union(a,c),d))]
new_row_index <- c(no_row, d, c, a)
#cols
a1 <- (1:100)[(nmtf_overlap$col_clusters[[1]][, 1])==1]
b1 <- (1:100)[(nmtf_overlap$col_clusters[[1]][, 2])==1]
c1 <- (1:100)[(nmtf_overlap$col_clusters[[1]][, 3])==1]
d1 <- (1:100)[(nmtf_overlap$col_clusters[[1]][, 4])==1]
no_col <-(1:100)[-sort(union(union(a1,b1),union(c1,d1)))]
new_col_index <- c(no_col, d1[!d1 %in% intersect(c1,d1)], 
        intersect(c1,d1), c1[!c1 %in% union(intersect(a1,c1), intersect(c1,d1))], 
        intersect(a1,c1), a1[!a1 %in% union(union(intersect(a1,b1), intersect(a1,d1)), intersect(a1,c1))], 
        intersect(a1,b1[!b1 %in% intersect(c1,b1)]), b1[!b1 %in% intersect(a1,b1)])
image((data_overlap$data_views[[1]])[new_row_index, new_col_index])

#view 2
x <- (1:200)[(nmtf_overlap$row_clusters[[2]][, 1])==1]
y <- (1:200)[(nmtf_overlap$row_clusters[[2]][, 2])==1]
z <- (1:200)[(nmtf_overlap$row_clusters[[2]][, 3])==1]
no_row <-(1:200)[-sort(union(union(x,y),z))]
new_row_index <- c(no_row, y, z, x)

x1 <- (1:250)[(nmtf_overlap$col_clusters[[2]][, 1])==1]
y1 <- (1:250)[(nmtf_overlap$col_clusters[[2]][, 2])==1]
z1 <- (1:250)[(nmtf_overlap$col_clusters[[2]][, 3])==1]
no_col <-(1:250)[-sort(union(union(x1,y1),z1))]
new_col_index <- c(no_col, y1[!y1 %in% intersect(y1,z1)], 
            intersect(y1,z1), z1[!z1 %in% union(intersect(y1,z1), intersect(x1,z1))], 
             intersect(x1,z1), x1[!x1 %in% union(intersect(x1,z1), intersect(x1,y1))])
new_col_index <- c(no_col, d1[!d1 %in% intersect(c1,d1)], 
        intersect(c1,d1), c1[!c1 %in% union(intersect(a1,c1), intersect(c1,d1))], 
        intersect(a1,c1), a1[!a1 %in% union(union(intersect(a1,b1), intersect(a1,d1)), intersect(a1,c1))], 
        intersect(a1,b1[!b1 %in% intersect(c1,b1)]), b1[!b1 %in% intersect(a1,b1)])

image((data_overlap$data_views[[2]])[new_row_index, new_col_index])


#using original relations
source("best_file_system/resNMTF_funcs.r")
nmtf_overlap2 <- restMultiNMTF_run(Xinput = data_overlap$data_views, 
            phi = phi_mat, relErr = "Sil")
source("best_file_system/evaluation_funcs.r")
nmtf_scores2 <- evaluate_simulation_comp(nmtf_overlap2$row_clusters , nmtf_overlap2$col_clusters,
 data_overlap$truth_rows, data_overlap$truth_cols, data_overlap$data_views, index = 2)


 #generate examples of biclustering 
 #now have overlapping biclusters
row_start <- list(c(1, 51, 151), c(1, 51, 151))
row_end <- list(c(75, 150, 200),  c(75, 150, 200))
col_start <- list(c(1, 31, 61),  c(1, 101, 151))
col_end <- list(c(40, 60, 100), c(120, 150, 250))
data_overlap <- multi_view(row_dims,
        col_dims, row_start, row_end, col_start, col_end, 
        noise=noise_level,row_same_shuffle=TRUE, col_same_shuffle=c(F, F), seed=10)
row_start <- list(c(1, 51), c(1, 51))
row_end <- list(c(75, 150),  c(75, 150))
col_start <- list(c(1, 31),  c(1, 101))
col_end <- list(c(40, 60), c(120, 150))
data_overlap <- multi_view(row_dims,
        col_dims, row_start, row_end, col_start, col_end, 
        noise=noise_level,row_same_shuffle=TRUE, col_same_shuffle=c(F, F), seed=10)
