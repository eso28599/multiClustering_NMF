#!/bin/bash
#PBS -N inc_phi4_res
#PBS -m a
#PBS -q medium
#PBS -o increasing_phi4/logs/results.out
#PBS -e increasing_phi4/logs/results.err

export sim_folder_name=increasing_phi4
export sim=phi

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla ${sim}_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt