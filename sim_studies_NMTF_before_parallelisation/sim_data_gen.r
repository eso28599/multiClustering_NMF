#file common to all sim studies
source("data_generation.r") #data gen files
source("sim_parameters.r") #parameters for simulation
#generate data and save in given path
save_data(n_repeats, row_cl_dims, col_cl_dims, path_to_data_folder, noise_level)
