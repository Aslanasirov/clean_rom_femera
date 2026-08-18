#!/bin/bash 
# Job Name
#$ -N FemeraSimulation
# Use current working directory
#$ -cwd
# Parallel Environment request.  Set your number of processors here
#$ -pe mpich 32
#$ -S /bin/bash
##### -l walltime=24:00:00
#$ -j y
#$ -o log.txt
ulimit -s unlimited















CompNodeNum=1
GrainNum=13
dGrainNum=1
CNum=32
nPhUCP=`echo $GrainNum / $CompNodeNum | bc`
icn=1
tmp=`echo $tmp*$nPhUCP | bc`
GrainID1=`echo $tmp+1 | bc`
GrainID2=`echo ${icn}*$nPhUCP | bc`
echo $GrainID1
echo $GrainID2
# Run Femera
./AutoRunFemera.sh $GrainNum $dGrainNum $GrainID1 $GrainID2 $CNum


