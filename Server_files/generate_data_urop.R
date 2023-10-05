library(mvtnorm)
library(Matrix)

generate_blocks_all <- function(row, column, rho=0.9, mean=5, mean_n=0, sd_n=1, overlap=FALSE, mul=FALSE, overlap_r=NULL, overlap_c=NULL, fullrank=TRUE, n_row_col=NULL, all_noise=FALSE) {
  
  set.seed(2)
  
  # For the purpose of generating matrix with all noise
  noise_mat <- list()
  
  if (all_noise==TRUE) {
    if (!is.null(n_row_col)) {
      for (i in length(n_row_col)) {
        noise_mat <- c(noise_mat, list(abs(matrix(rnorm(n = n_row_col[[i]][[1]] * n_row_col[[i]][[2]], mean = mean_n, sd = sd_n), nrow = n_row_col[[i]][[1]], ncol = n_row_col[[i]][[2]]))))
      }
    }
    else {
      for (i in 1:3) {
        noise_mat <- c(noise_mat, list(abs(matrix(rnorm(n = 22500, mean = mean_n, sd = sd_n), nrow = 150, ncol = 150))))
      }
    }
    return(noise_mat)
  }
  
  # Test the row and column
  if (length(row) != length(column) || any(sapply(seq_along(row), function(i) length(row[[i]]) != length(column[[i]])))) {
    stop("Error: Lengths of row and column must be the same and the lengths of corresponding elements must be the same.")
  }
  
  result = list()
  shuffle = list()
  indices = list()
  biclust_list = list()
  row_c_list = list()
  col_c_list = list()
  
  for (i in 1:length(row)) {
    
    result_no_noise <- NULL
    result_no_noise_m <- matrix(0, nrow = sum(unlist(row[[i]])), ncol = sum(unlist(column[[i]])))
    ol_mat <- NULL
    
    ol_r <- 0
    ol_c <- 0
    
    nol_r <- 0
    nol_c <- 0
    
    row_c <- matrix(0, nrow = sum(unlist(row[[i]])), ncol = length(row[[i]]))
    col_c <- matrix(0, nrow = sum(unlist(column[[i]])), ncol = length(row[[i]]))
    r_index <- 0
    c_index <- 0
    
    for (j in 1:length(row[[i]])) {
      r <- row[[i]][[j]]
      c <- column[[i]][[j]]
      H <- abs(outer(1:c, 1:c, "-"))
      sigma <- rho^H
      block <- rmvnorm(r, mean=rep(mean, c), sigma=sigma)
      
      # Original block structure
      if (overlap==FALSE) {
        if (j == 1) {
          result_no_noise <- block
          row_c[1:r, 1] <- 1
          col_c[1:c, 1] <- 1
          r_index <- r
          c_index <- c
        }
        else{
          result_no_noise <- bdiag(result_no_noise, block)
          row_c[(r_index+1):(r_index+r), j] <- 1
          col_c[(c_index+1):(c_index+c), j] <- 1
          r_index <- r_index+r
          c_index <- c_index+c
        }
      }
      # Overlapping structure
      else {
        
        if (j == 1) {
          result_no_noise_m[1:r, 1:c] <- block
          row_c[1:r, 1] <- 1
          col_c[1:c, 1] <- 1
          r_index <- r
          c_index <- c
        }
        
        else {
          result_no_noise_m[(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+1):(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+r), (sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+1):(sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+c)] <- block
          if (mul) {
            result_no_noise_m[(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+1):(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+tail(nol_r, 1)), (sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+1): (sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+tail(nol_c, 1))] <- ol_mat*block[1:nrow(ol_mat), 1:ncol(ol_mat)]
          }
          else {
            result_no_noise_m[(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+1):(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+tail(nol_r, 1)), (sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+1): (sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+tail(nol_c, 1))] <- ol_mat+block[1:nrow(ol_mat), 1:ncol(ol_mat)]
          }
          
          # Generate the biclusters
          row_c[(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+1):(sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+r), j] <- 1
          col_c[(sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+1):(sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+c), j] <- 1
          r_index <- (sum(unlist(row[[i]])[1:(j-1)])-sum(nol_r)+r)
          c_index <- (sum(unlist(column[[i]])[1:(j-1)])-sum(nol_c)+c)
        }
        if (j==length(row[[i]])) {
          break
        }
        
        r_m <- r-tail(nol_r, 1)
        c_m <- c-tail(nol_c, 1)
        if (!is.null(overlap_r)) {
          ol_r <- overlap_r[[j]]
        }
        else {
          n <- min(r_m, row[[i]][[j+1]])
          ol_r <- sample.int(n-1, 1, prob = 1+0.1*(n-(1:(n-1))))+1
        }
        if (!is.null(overlap_c)) {
          ol_c <- overlap_c[[j]]
        }
        else {
          m <- min(c_m, column[[i]][[j+1]])
          ol_c <- sample.int(m-1, 1, prob = 1+0.1*(m-(1:(m-1))))+1
        }
        ol_mat <- block[(r-ol_r+1):r, (c-ol_c+1):c]
        nol_r <- c(nol_r, ol_r)
        nol_c <- c(nol_c, ol_c)
      }
    }
    
    
    if (!identical(result_no_noise_m, matrix(0, nrow = sum(unlist(row[[i]])), ncol = sum(unlist(column[[i]]))))) {
      result_no_noise <- result_no_noise_m
    }
    
    
    # Get rid of all empty rows and columns
    empty_rows <- rowSums(result_no_noise) == 0
    empty_cols <- colSums(result_no_noise) == 0
    result_no_noise <- result_no_noise[!empty_rows, !empty_cols]
    
    row_c <- row_c[1:r_index, ]
    col_c <- col_c[1:c_index, ]
    
    # Add empty rows and columns if is needed
    if (fullrank==FALSE) {
      if (!is.null(n_row_col)) {
        result_no_noise <- rbind(as.matrix(result_no_noise), matrix(0, nrow = n_row_col[1]-nrow(result_no_noise), ncol = ncol(result_no_noise)))
        result_no_noise <- cbind(as.matrix(result_no_noise), matrix(0, nrow = nrow(result_no_noise), ncol = n_row_col[2]-ncol(result_no_noise)))
        row_c <- rbind(row_c, matrix(0, nrow = n_row_col[1]-nrow(row_c), ncol = ncol(row_c)))
        col_c <- rbind(col_c, matrix(0, nrow = n_row_col[2]-ncol(col_c), ncol = ncol(col_c)))
      }
      else {
        result_no_noise <- rbind(as.matrix(result_no_noise), matrix(0, nrow = sample.int(20, 1), ncol = ncol(result_no_noise)))
        result_no_noise <- cbind(as.matrix(result_no_noise), matrix(0, nrow = nrow(result_no_noise), ncol = sample.int(20, 1)))
        row_c <- rbind(row_c, matrix(0, nrow = sample.int(20, 1), ncol = ncol(row_c)))
        col_c <- rbind(col_c, matrix(0, nrow = sample.int(20, 1), ncol = ncol(col_c)))
      }
    }
    
    result_no_noise <- as.matrix(result_no_noise)
    
    noise <- matrix(rnorm(n = nrow(result_no_noise) * ncol(result_no_noise), mean = mean_n, sd = sd_n), nrow = nrow(result_no_noise), ncol = ncol(result_no_noise))
    new_result <- abs(result_no_noise+noise)
    
    result <- c(result, list(new_result))
    
    
    n_rows <- nrow(new_result)
    n_cols <- ncol(new_result)
    
    shuffled_matrix <- new_result
    
    # Randomly select two distinct indices for rows or columns
    indices_r <- sample(n_rows, n_rows, replace = FALSE)
    # Swap the selected rows
    shuffled_matrix <- shuffled_matrix[indices_r, ]
    
    # For columns 
    indices_c <- sample(n_cols, n_cols, replace = FALSE)
    shuffled_matrix <- shuffled_matrix[, indices_c]
    
    
    shuffle <- c(shuffle, list(shuffled_matrix))
    indices <- c(indices, list(indices_r, indices_c))
    row_c_list <- c(row_c_list, list(row_c))
    col_c_list <- c(col_c_list, list(col_c))
  }
  
  
  final <- list("Original" = result, "Shuffled" = shuffle, "Indices" = indices, "row_c" = row_c_list, "col_c" = col_c_list)
  return(final)
}
```
