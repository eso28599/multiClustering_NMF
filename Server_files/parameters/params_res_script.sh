#!/bin/bash
#PBS -N R_job
#PBS -m a
#PBS -q jumbo
#PBS -o NMTF_params/logs/results.out
#PBS -e NMTF_params/logs/results.err

export sim_folder_name=NMTF_params

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla params_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt