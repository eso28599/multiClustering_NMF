#!/bin/bash
#PBS -N inc_views_res
#PBS -m a
#PBS -q medium
#PBS -o increasing_views/logs/results.out
#PBS -e increasing_views/logs/results.err

export sim_folder_name=increasing_views

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla all_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt