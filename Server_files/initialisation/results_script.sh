#!/bin/bash
#PBS -N R_job
#PBS -m a
#PBS -q jumbo
#PBS -o NMTF_initialisation/logs/results.out
#PBS -e NMTF_initialisation/logs/results.err

export sim_folder_name=NMTF_initialisation

cd ${PBS_O_WORKDIR}/${sim_folder_name}
ls data/ > repeat_folders.txt
#move back into original folder
cd ${PBS_O_WORKDIR}


#compile results
Rscript --vanilla init_results.r  ${sim_folder_name}  ${sim_folder_name}/repeat_folders.txt