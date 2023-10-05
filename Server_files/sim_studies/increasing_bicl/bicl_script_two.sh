#!/bin/bash
#PBS -N inc_bicl_res
#PBS -m a
#PBS -q medium
#PBS -o increasing_bicl/logs/results.out
#PBS -e increasing_bicl/logs/results.err

export sim_folder_name=increasing_bicl

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla all_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt
