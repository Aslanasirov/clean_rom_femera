#!/bin/bash

# Job Name
#$ -N P5e6
# Use current working directory
#$ -cwd
###$ -V
# Parallel Environment request.  Set your number of processors here
#$ -pe mpich 1

# Run job through bash shell
#$ -S /bin/bash
# INPUT: RUN INTEL ONEAPI SCRIPT TO SET UP FORTRAN PATH
source /home/anasirov/intel/oneapi/setvars.sh intel64 --force



#=============================================================

time1="$(date -u +%s.%N)"
# EXECUTE RUNSPARSE CODE
./RunSparse.sh >> ./runsparselog

time2="$(date -u +%s.%N)"
elapsed="$(bc <<<"$time2-$time1")"
echo "***** Execution time of RunSparse is "
echo $elapsed

