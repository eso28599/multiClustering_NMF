#!/bin/bash
#PBS -N inc_phi3_res
#PBS -m a
#PBS -q medium
#PBS -o increasing_phi3/logs/results.out
#PBS -e increasing_phi3/logs/results.err

export sim_folder_name=increasing_phi3
export sim=phi

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla ${sim}_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt