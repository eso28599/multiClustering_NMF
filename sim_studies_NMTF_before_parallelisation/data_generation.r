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

multi_view <- function(rowClusters,colClusters,noise=10,seed=FALSE){
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

    #generate data and store
    for(i in 1:n_views){
        data_i <- one_view(rowClusters[[i]],colClusters[[i]],noise)
        X_trial[[i]] <- data_i$view
        true_row_clusterings[[i]] <- data_i$truth_row
        true_col_clusterings[[i]] <- data_i$truth_col
    }

    return(list(data_views= X_trial, truth_rows = true_row_clusterings, truth_cols = true_col_clusterings))
}

save_data <- function(N_repeats,rowClusters,colClusters,file_path, noise=10,seed=FALSE){
    for (i in 1:N_repeats){
        #can change the noise parameter here for level of noise in views
        our_data <- multi_view(rowClusters,colClusters,noise)
        #save data as a file in current directory 
        #export each data frame to separate sheets in same Excel file
        data_name <- paste0(paste0(paste0(file_path,'/data/repeat'),i),'.xlsx')
        row_clusts_name <- paste0(paste0(paste0(file_path,'/true_rows/repeat'),i),'.xlsx')
        col_clusts_name <- paste0(paste0(paste0(file_path,'/true_cols/repeat'),i),'.xlsx')
        openxlsx::write.xlsx(our_data$data_views, file = data_name)
        openxlsx::write.xlsx(our_data$truth_rows, file = row_clusts_name)
        openxlsx::write.xlsx(our_data$truth_rows, file = col_clusts_name)
    }
}
