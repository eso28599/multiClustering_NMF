# common functions 
library(rio)
# generate names of files 
generate_file_names <- function(n_repeats,path_to_data_folder,file_sub_path){
    filenames <- sapply(1:n_repeats, function(x) {
                    paste0(paste0(paste0(path_to_data_folder,file_sub_path),paste0("/repeat",i)),'.xlsx')
                })
return(filenames)
}
#import as list of matrices instead of as a list
import_matrix <- function(filename){
    return(lapply(import_list(filename), function(x) as.matrix(x)))
}