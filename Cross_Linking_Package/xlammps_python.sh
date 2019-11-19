#! /bin/bash
# 
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M hnalmahm@mtu.edu
#$ -m abes
#$ -pe mpich_unstaged 16
#$ -q long.q
#$ -hard -l mem_free=2G
#$ -hard -l lammps_lic=.06250000
#$ -notify

# Necessary variables
source /share/apps/bin/bashrc
module list
module purge
module load intel/2016.1
module load lammps/2016.02.29-CPU
module list
which python

# Input/Output files
export INPUT_FOLDER="${PWD}"
export INPUT_FILE="xlinkCN.script"
export OUTPUT_FILE="outCNxlink.out"
export ARRAY_JOB=""

# Run LAMMPS (parallel) + Python post-processing
export IMAX="2"

# Required crosslink percentage (max xlink %)
export xlmax="80"     

###############################
######## LOOP START ###########
###############################

i=1
while [ ${i} -le ${IMAX} ]
do

  ####################
  ###### LAMMPS ######
  ####################

  cd ${INPUT_FOLDER}/lammps
  mpirun -n ${NSLOTS} -machine ${TMP}/machines ${LAMMPS}/src/lmp_mtu_parallel -log ${INPUT_FOLDER}/lammps/${OUTPUT_FILE}${ARRAY_JOB} -in ${INPUT_FOLDER}/lammps/${INPUT_FILE}${ARRAY_JOB}
  cp data.xlinkNC ${INPUT_FOLDER}/python/

  ##################
  ##### PYTHON #####
  ##################

  cd ${INPUT_FOLDER}/python
  python xupdate.py
  cp x_update.mol ${INPUT_FOLDER}/lammps/
  

  ###################################
  ### Required Crosslink  Checker ###
  ###################################
 
  xlper=$(grep 'CROSSlinkpercent' `ls -dt ${INPUT_FOLDER}/xlammps_python.sh.o* | head -1` | tail -1 | awk -F' ' '{print $2}' | awk -F'.' '{print $1}')
  if [ ${xlper} -ge ${xlmax} ]
  then
   break  # Skip entire rest of loop.
  fi

  # Increase i by one
  i=$(expr ${i} + 1)
done

# Unload modules
module unload lammps/2016.02.29-CPU
module unload intel/2016.1
module list
