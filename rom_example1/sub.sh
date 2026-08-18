#!/bin/bash
# THIS SCRIPT RUNS ThE ROM ON THE CLUSTER
# CALLED FROM ./caller.sh from Simulation folder

cwd=$(pwd)
echo pwd: $cwd

rm runsparselog
rm simlog

qsub -cwd -j y -o simlog -l h=node013 inter.sh


