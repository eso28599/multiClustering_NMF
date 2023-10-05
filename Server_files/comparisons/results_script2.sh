#!/bin/bash
#PBS -N R_job
#PBS -m a
#PBS -q medium
#PBS -o scen_2/logs/results.out
#PBS -e scen_2/logs/results.err

export sim_folder_name=scen_2

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla all_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt