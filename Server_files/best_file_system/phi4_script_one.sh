#!/bin/bash
#PBS -N increasing_phi4
#PBS -m a
#PBS -q medium
#PBS -t 1-50
#PBS -o increasing_phi4/logs/test_job.out
#PBS -e increasing_phi4/logs/test_job.err

export sim_folder_name=increasing_phi4
export sim=phi
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
Rscript --vanilla ${sim}_data_gen.r  ${sim_folder_name} $I

#now analyse in R
Rscript --vanilla ${sim}_methods_r.r  ${sim_folder_name} $I

#analyse in python
#python3 methods_python.py ${sim_folder_name} $I

#evaluate results in R
Rscript --vanilla ${sim}_eval.r  ${sim_folder_name} $I
