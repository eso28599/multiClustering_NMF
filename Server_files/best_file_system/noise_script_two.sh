#!/bin/bash
#PBS -N inc_noise
#PBS -m a
#PBS -q medium
#PBS -o increasing_noise/logs/results.out
#PBS -e increasing_noise/logs/results.err

export sim_folder_name=increasing_noise
export sim=noise

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla ${noise}_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt
