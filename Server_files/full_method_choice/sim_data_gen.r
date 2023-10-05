args = commandArgs(trailingOnly = TRUE)
path_to_sim_folder = as.character(args[1])
batch_folder = as.character(args[2])
#make the simulation folder name the argument
#create arg which is something like NMTF_initialisation
#file common to all sim studies
source(paste0(path_to_sim_folder, "/sim_parameters.r")) #parameters for simulation
source("data_generation.r") #data gen files
#generate data and save in given path
path_to_data_folder <- paste0(paste0(path_to_sim_folder, "/data/"),batch_folder)
save_data(row_cl_dims, col_cl_dims, path_to_data_folder, row_same_shuffle, col_same_shuffle, noise_level)
