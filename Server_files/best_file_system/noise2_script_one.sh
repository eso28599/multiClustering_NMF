#!/bin/bash
#PBS -N increasing_noise2
#PBS -m a
#PBS -q medium
#PBS -t 1-50
#PBS -o increasing_noise2/logs/test_job.out
#PBS -e increasing_noise2/logs/test_job.err

export sim_folder_name=increasing_noise2
export sim=noise
export i=${PBS_ARRAYID}
export I=`echo $i | awk '{printf "%3.3d", $1}'`


cd ${PBS_O_WORKDIR}/${sim_folder_name}/data
#mkdir $I
if [ ! -d "$I" ]; then
  mkdir $I
  cd $I
  for j in 1 10 20 30 40 50 60 70 80 90 100
  do
    mkdir res_nmtf_$j
  done
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
