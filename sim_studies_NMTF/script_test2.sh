#!/bin/bash
#PBS -N R_job
#PBS -m e
#PBS -q medium
#PBS -o test_job.out
#PBS -e test_job.err

export sim_folder_name=NMTF_initialisation
export i=$SGE_TASK_ID
export I=`echo $i | awk '{printf "%3.3d", $1}'`
#cd ${HOME}/data_integration/sim_studies_NMTF

cd ${PBS_O_WORKDIR}/${sim_folder_name}/data
mkdir $I
#move back into original folder
cd ${PBS_O_WORKDIR}

#generate data
Rscript --vanilla sim_data_gen.r  ${sim_folder_name} $I

#now analyse in R
Rscript --vanilla init_standard_kmeans.r  ${sim_folder_name} $I
