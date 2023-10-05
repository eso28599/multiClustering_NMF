#!/bin/bash
#PBS -N R_job
#PBS -m a
#PBS -q medium
#PBS -t 1-1000
#PBS -o k_end_choice/logs/test_job.out
#PBS -e k_end_choice/logs/test_job.err

export sim_folder_name=k_end_choice
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/${sim_folder_name}/data
mkdir $I
#move back into original folder
cd ${PBS_O_WORKDIR}

#generate data
Rscript --vanilla sim_data_gen.r  ${sim_folder_name} $I

#now analyse in R
Rscript --vanilla init_standard_kmeans.r  ${sim_folder_name} $I

#evaluate results in R
Rscript --vanilla init_eval.r  ${sim_folder_name} $I