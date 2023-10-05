library(openxlsx)
library(MASS)
library(Matrix)
##### synthetic data generation

#preparation of linked dataviews
# function to generate a dataview for specified bicluster dimensions
# same number of row and column clusters 
# row clusters remain common across all views, columns change 

 one_view_adv <- function(row_dims, col_dims, row_start, row_end,  col_start, col_end, overlap = FALSE, noise=10){
    #'row_dims: vector of sizes of each row cluster in this view
    #'col_dims: vector of sizes of each row cluster in this view
    #'noise: variance of the noise added to views
    #' 
    #'view: matrix of this dataview
    #'truth_row: vector indicating membership of row clusters
    #'truth_col: vector indicating membership of col clusters

    # generate noise
    X_noise <- mvrnorm(n = row_dims, mu = rep(0, col_dims), Sigma = diag(noise, col_dims))

    
    # Create a list as input for block-diagonal data generation
    # list length of no of clusters in this first view
    no_row_cl <- length(row_start)
    no_col_cl <- length(col_start)

    # generate a mvn with n=number of individuals of view and no of features
    # equal to number of features of the view
    # mean of 5 for each column
    # covariance matrix is identity - each feature is independent
    #define vectors indicating true row/column cluster membership for each view
    true_row <- matrix(0, nrow = row_dims, ncol = no_row_cl)
    true_col <- matrix(0, nrow = col_dims, ncol = no_col_cl)
    X_view <- matrix(0, nrow = row_dims, ncol = col_dims)
    for (i in 1:no_row_cl){
        n_r <- (row_end[[i]] - row_start[[i]] + 1)
        n_c <- (col_end[[i]] - col_start[[i]] + 1)
        X_view[(row_start[i]):(row_end[[i]]), (col_start[[i]]):(col_end[[i]])] <- X_view[(row_start[[i]]):(row_end[[i]]),
                 (col_start[[i]]):(col_end[[i]])] + 
                    mvrnorm(n = n_r, mu = rep(5, n_c), Sigma = diag(n_c))
        true_row[(row_start[[i]]):(row_end[[i]]), i] <- 1
        true_col[(col_start[[i]]):(col_end[[i]]), i] <- 1
    }
    # add noise to first view
    X_view <- abs(X_view) + abs(X_noise)
    return(list(view = X_view, truth_row = true_row, truth_col = true_col))
}

multi_view <- function(row_dims, col_dims, row_start, row_end,  col_start, col_end, noise=10,row_same_shuffle=TRUE, col_same_shuffle=TRUE, seed=FALSE){
    #'rowClusters: n length list of vectors of row cluster sizes in each view
    #'rowClusters: n length list of vectors of column cluster sizes in each view
    #'seed: logical indicator, default is FALSE, if true sets seed so same data is generated each time
    #' 
    #'data_views: n length list of data views
    #'truth_rows: n length list of vectors of true row cluster membership for each view
    #'col_rows: n length list of vectors of true column cluster membership for each view
    
    if(seed){
        set.seed(20)
    }
    n_views <- length(row_dims)
    # Introduce simulated data-views- where we store the views
    X_trial <- vector("list", length = n_views)
    #list to store row clusterings for each dataview
    true_row_clusterings <- vector("list", length = n_views)
    true_col_clusterings <- vector("list", length = n_views) 
    #initialise index to shuffle rows and columns 
    #new row/col index
    #if shuffle is true - dims of each view must be the same and represent same objects
    if(row_same_shuffle){
        new_row_ind <- sample(row_dims[1])
    }
    if(col_same_shuffle){
        new_col_ind <- sample(col_dims[1])
    }
    #generate data and store
    for(i in 1:n_views){
        if(!row_same_shuffle){
            new_row_ind <- sample(row_dims[i])
        }
        if(!col_same_shuffle){
            new_col_ind <- sample(col_dims[i])
        }
        data_i <- one_view_adv(row_dims[i], col_dims[i],
                 row_start[[i]], row_end[[i]],  col_start[[i]],
                  col_end[[i]], noise)
        X_trial[[i]] <- (data_i$view)[new_row_ind, new_col_ind]
        true_row_clusterings[[i]] <- (data_i$truth_row)[new_row_ind, ]
        true_col_clusterings[[i]] <- (data_i$truth_col)[new_col_ind, ]
    }

    return(list(data_views= X_trial, truth_rows = true_row_clusterings, truth_cols = true_col_clusterings))
}

save_data <- function(row_dims, col_dims, row_start, row_end,  col_start, col_end, file_path,row_same_shuffle=TRUE,col_same_shuffle=TRUE, noise = 10){
        #can change the noise parameter here for level of noise in views
        data <- multi_view(row_dims, col_dims, row_start, row_end,  col_start, col_end, noise, row_same_shuffle , col_same_shuffle )
        #save data as a file in given directory
        #export each data frame to separate sheets in same Excel file
        openxlsx::write.xlsx(data$data_views, file = paste0(file_path, "/data.xlsx")) # nolint
        openxlsx::write.xlsx(data$truth_rows, file = paste0(file_path, "/true_rows.xlsx")) # nolint: line_length_linter.
        openxlsx::write.xlsx(data$truth_cols, file = paste0(file_path, "/true_cols.xlsx")) # nolint: line_length_linter.
}
