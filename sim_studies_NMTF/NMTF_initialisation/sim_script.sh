#!/bin/bash
#PBS -N standard_submission_script
#PBS -m be
#PBS -q standard
#PBS -t 1-2
#PBS -o ${HOME}/data_integration/sim_studies_NMTF/logs
#PBS -e ${HOME}/data_integration/sim_studies_NMTF/logs

export sim_folder_name = "NMTF_initialisation"
export i= $SGE_TASK_ID
export I=`echo $i | awk '{printf "%3.3d", $1}'`

cd ${HOME}/data_integration/sim_studies_NMTF/${sim_folder_name}

#Script for carrying out simulation - save as name of simulation
#2. Create empty folder with the name 'sim_name_of_sim_data' following subfolders 


cd data #move into the data folder
mkdir $I #make new folder for this repeat
cd ${HOME}/data_integration/sim_studies_NMTF #move back to where the files are for this repeat
#2. Generate data and store in current folder


/usr/bin/R --vanilla < sim_data_gen.r > $sim_folder_name  $I output.$I


#3. Run analysis in R, saving clusters
#4. Run evaluation of results and save 

