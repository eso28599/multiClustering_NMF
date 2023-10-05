library(openxlsx)
library(MASS)
library(Matrix)
##### synthetic data generation

#preparation of linked dataviews
# function to generate a dataview for specified bicluster dimensions
# same number of row and column clusters 
# row clusters remain common across all views, columns change 

one_view <- function(row_dims,col_dims,noise=10){
    #'row_dims: vector of sizes of each row cluster in this view
    #'col_dims: vector of sizes of each row cluster in this view
    #'noise: variance of the noise added to views
    #' 
    #'view: matrix of this dataview
    #'truth_row: vector indicating membership of row clusters
    #'truth_col: vector indicating membership of col clusters
    n_row <- length(row_dims)
    n_col <- length(col_dims)


    # Create a list as input for block-diagonal data generation
    # list length of no of clusters in this first view
    inputList <- vector("list", length = max(n_row, n_col))
    # generate a mvn with n=number of individuals of view and no of features
    # equal to number of features of the view
    # mean of 5 for each column
    # covariance matrix is identity - each feature is independent 
    #define vectors indicating true row/column cluster membership for each view
    true_row <- c()
    true_col <- c()
    for(i in 1: n_row ){
        inputList[[i]] <- mvrnorm(n =row_dims[i], mu = rep(5,col_dims[i]), Sigma = diag(col_dims[i]))
        true_row <- c(true_row,rep(i,row_dims[i]))
        true_col <- c(true_col,rep(i,col_dims[i]))
    }

    X_view <- as.matrix(bdiag(inputList)) 
    # generate noise, MVN with no correlation, zero mean and variance of 10
    X_noise <- mvrnorm(n = nrow(X_view), mu = rep(0, ncol(X_view)), Sigma = diag(noise,ncol(X_view)))
    # add noise to first view
    X_view <- X_view + X_noise 
    return(list(view= X_view, truth_row = true_row, truth_col = true_col))

}

multi_view <- function(rowClusters,colClusters,noise=10,row_same_shuffle=TRUE, col_same_shuffle=TRUE, seed=FALSE){
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
    n_views <- length(rowClusters)
    # Introduce simulated data-views- where we store the views
    X_trial <- vector("list", length = n_views)
    #list to store row clusterings for each dataview
    true_row_clusterings <- vector("list", length = n_views) 
    true_col_clusterings <- vector("list", length = n_views) 
    #initialise index to shuffle rows and columns 
    #new row/col index
    #if shuffle is true - dims of each view must be the same and represent same objects
    if(row_same_shuffle){
        new_row_ind <- sample(sum(rowClusters[[1]]))
    }
    if(col_same_shuffle){
        new_col_ind <- sample(sum(colClusters[[1]]))
    }
    #generate data and store
    for(i in 1:n_views){
        if(!row_same_shuffle){
            new_row_ind <- sample(sum(rowClusters[[i]]))
        }
        if(!col_same_shuffle){
            new_col_ind <- sample(sum(colClusters[[i]]))
        }
        data_i <- one_view(rowClusters[[i]],colClusters[[i]],noise)
        X_trial[[i]] <- (data_i$view)[new_row_ind, new_col_ind]
        true_row_clusterings[[i]] <- (data_i$truth_row)[new_row_ind]
        true_col_clusterings[[i]] <- (data_i$truth_col)[new_col_ind]
    }

    return(list(data_views= X_trial, truth_rows = true_row_clusterings, truth_cols = true_col_clusterings))
}

save_data <- function(row_clusters, col_clusters, file_path,row_same_shuffle=TRUE,col_same_shuffle=TRUE, noise = 10){
        #can change the noise parameter here for level of noise in views
        data <- multi_view(row_clusters, col_clusters, noise, row_same_shuffle , col_same_shuffle )
        #save data as a file in given directory
        #export each data frame to separate sheets in same Excel file
        openxlsx::write.xlsx(data$data_views, file = paste0(file_path, "/data.xlsx")) # nolint
        openxlsx::write.xlsx(data$truth_rows, file = paste0(file_path, "/true_rows.xlsx")) # nolint: line_length_linter.
        openxlsx::write.xlsx(data$truth_cols, file = paste0(file_path, "/true_cols.xlsx")) # nolint: line_length_linter.
}
