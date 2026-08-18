#!/bin/bash

# SET UP PYTHON PATH
python3=/home/anasirov/mylibs/python3/bin/python3

echo Got $NSLOTS processors.
echo Machines:
echo JOBNAME: $JOBNAME
echo JOBN0: $JOBN0
echo JOB_ID: ${JOB_ID}
cat $TMPDIR/machines

# get the current working directory
cwd=$(pwd)
echo pwd: $cwd

sleep 0.0

echo "============="
echo " MODIFY SIM FILES TO HAVE CORRECT INPUTS"
echo "============="
GrainNum=13 # INPUT: NUMBER OF GRAINS
$python3 ModifySimFiles.py $GrainNum

# Define the Abaqus particulars
JOB_NAME=tension # INPUT: SPECIFY ABAQUS INPUT FILE NAME
INPUT_FILE=${cwd}/Simulation/${JOB_NAME}.inp

# SPECIFY ABAQUS PATH
ABAQUS_ARGS="user=${cwd}/Simulation/UMAT.f -inter"
ABAQUS="/cm/shared/apps/SIMULIA/var/Commands/abaqus"

## Run the job
cd $cwd/Simulation
rm ${JOB_NAME}.lck
ulimit -s unlimited
$ABAQUS job=${JOB_NAME} input=${INPUT_FILE} ${ABAQUS_ARGS}
