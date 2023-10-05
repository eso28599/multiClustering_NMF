#!/bin/bash
#PBS -N R_job
#PBS -m a
#PBS -q medium
#PBS -t 1-100
#PBS -o NMTF_params_diff_scen_s/logs/test_job.out
#PBS -e NMTF_params_diff_scen_s/logs/test_job.err

export sim_folder_name=NMTF_params_diff_scen_s
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/${sim_folder_name}/data
mkdir $I
#move back into original folder
cd ${PBS_O_WORKDIR}

#generate data
Rscript --vanilla sim_data_gen.r  ${sim_folder_name} $I

#now analyse in R
Rscript --vanilla params_methods.r  ${sim_folder_name} $I

#evaluate results in R
Rscript --vanilla params_eval.r  ${sim_folder_name} $I
