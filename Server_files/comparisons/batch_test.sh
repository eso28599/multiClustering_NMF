#!/bin/bash
#PBS -N increasing_views
#PBS -m a
#PBS -q medium
#PBS -t 1-25
#PBS -o scen_1/logs/test_job.out
#PBS -e scen_1/logs/test_job.err

export sim_folder_name=increasing_views
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/${sim_folder_name}/data
#mkdir $I
if [ ! -d "$I" ]; then
  mkdir $I
fi
#move back into original folder
cd ${PBS_O_WORKDIR}

#generate data
Rscript --vanilla sim_data_gen.r  ${sim_folder_name} $I

#now analyse in R
Rscript --vanilla methods_r.r  ${sim_folder_name} $I

#analyse in python
#python3 methods_python.py ${sim_folder_name} $I

#evaluate results in R
Rscript --vanilla all_eval.r  ${sim_folder_name} $I
