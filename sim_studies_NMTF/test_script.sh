#!/bin/sh
#PBS -N standard_submission_script
#PBS -m be
#PBS -q standard
#PBS -o /logs
#PBS -e /logs


#  test_script.sh
#  
#
#  Created by Ella Orme on 24/03/2023.
#  

#export sim_folder_name = "NMTF_initialisation"
#export i= $SGE_TASK_ID
#export I=`echo $i | awk '{printf "%3.3d", $1}'`

export I="1"

#cd ${PBS_O_WORKDIR}/NMTF_initialisation

#Script for carrying out simulation - save as name of simulation
#2. Create empty folder with the name 'sim_name_of_sim_data' following subfolders


cd ${PBS_O_WORKDIR}/data #move into the data folder
mkdir $I #make new folder for this repeat
#cd ${HOME}/data_integration/sim_studies_NMTF #move back to where the files are for this repeat
#2. Generate data and store in current folder

cd ${PBS_O_WORKDIR}
#/usr/bin/R --vanilla < sim_data_gen.r > $sim_folder_name  $I output.$I
/usr/bin/R --vanilla < sim_data_gen.r > /NMTF_initialisation  $I output.$I
