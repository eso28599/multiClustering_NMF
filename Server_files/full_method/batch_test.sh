#!/bin/bash
#PBS -N R_job
#PBS -m a
#PBS -q jumbo
#PBS -t 1-25
#PBS -o NMTF_initialisation/logs/test_job.out
#PBS -e NMTF_initialisation/logs/test_job.err

export sim_folder_name=NMTF_initialisation
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