#!/bin/bash 
#$ -o amp.txt
#$ -N AssembleAMP
#$ -j y
#$ -cwd
# Parallel Environment request.  Set your number of processors here
#$ -pe mpich 1
#$ -l h=node013
#$ -S /bin/bash
# DEFINE INPUTS
# Define parameters
CompNodeNum=1
GrainNum=13
dGrainNum=1

 
$HOME/mylibs/python3/bin/python3 AssembleAMP_modified.py $GrainNum $dGrainNum $CompNodeNum

$HOME/mylibs/python3/bin/python3 NewAMPCalStrain1Layer.py $GrainNum
