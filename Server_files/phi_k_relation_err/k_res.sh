#!/bin/bash
#PBS -N R_job
#PBS -m a
#PBS -q jumbo
#PBS -o k_phi/logs/results.out
#PBS -e k_phi/logs/results.err

export sim_folder_name=k_phi

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla params_results_sd.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt