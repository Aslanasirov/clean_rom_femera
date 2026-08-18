#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N prepare
#$ -j y
#$ -o log.txt
#$ -l h=node013
# ====== This file is employed to test the parallelization capability of Femera ========

# get the current working directory
cwd=$(pwd)
echo pwd: $cwd

# Define parameters
CompNodeNum=1 # Input: Number of computational nodes
GrainNum=13 # Input: number of grains
dGrainNum=1 # Input: just set to 1
R=1 #INPUT:  VOXEL SIZE (RESOLUTION)
Prefix='example1' #INPUT: PREFIX

# CUBIC ELASTIC CONSTANTS
# ALUMINUM 
#C11bcc=108200.0
#C12bcc=61300.0
#C44bcc=28500.0
# INCONEL
C11bcc=243300.0 
C12bcc=156700.0
C44bcc=117800.0


L1=10.0; L2=10.0; L3=10.0; # INPUT: SIZE OF THE RVE
CNum=32 # INPUT: NUMBER OF OPENMP CORES TO BE USED BY FEMERA

# python path
python3=$HOME/mylibs/python3/bin/python3


#==========================================================================
c1=$C11bcc ; c2=$C11bcc ; c3=$C11bcc ; 
c4=$C12bcc ; c5=$C12bcc ; c6=$C12bcc ; 
c7=$C44bcc ; c8=$C44bcc ; c9=$C44bcc ; 


# ======================== GENERATE MESH FILES ====================================
[ -d "./Partitions" ] && rm -r Partitions
mkdir Partitions
mkdir Partitions/GrainEls
[ -d "./PBC" ] && rm -r PBC
mkdir PBC


$python3 ./PartitionSplit.py $GrainNum $Prefix
$python3 ./PartitionPBCs.py $GrainNum $L1 $L2 $L3 $R $Prefix
# ======================== Convert mesh files into .fmr files ===========================
CNumber=1 #
[ -d "./femera-Neper" ] && rm -r femera-Neper
scp -r ./femera-Backup ./femera-Neper
cd ./femera-Neper
# Modify RunNeper.sh file
scp ../RunNeper.sh RunNeper.sh
./RunNeper.sh $GrainNum $dGrainNum $CNumber $c1 $c2 $c3 $c4 $c5 $c6 $c7 $c8 $c9 $L1 $L2 $L3

$python3 ../FixBCnElast.py $GrainNum $dGrainNum $c1 $c2 $c3 $c4 $c5 $c6 $c7 $c8 $c9 $L1 $L2 $L3
cd ..

# ======================== Modify Run.sh ===========================
sed -i "7s/.*/#$ -pe mpich $CNum/" ./Run.sh
sed -i "28s/.*/CompNodeNum=$CompNodeNum/" ./Run.sh
sed -i "29s/.*/GrainNum=$GrainNum/" ./Run.sh
sed -i "30s/.*/dGrainNum=$dGrainNum/" ./Run.sh
sed -i "31s/.*/CNum=$CNum/" ./Run.sh

# ============= Create folder and subfolders =============
SubmitFolder='Jobs'
[ -d "./$SubmitFolder" ] && rm -r $SubmitFolder

# Modify Submit.sh
sed -i "3s/.*/CompNodeNum=$CompNodeNum/" ./Submit.sh

# Prepare files for each computational node
mkdir $SubmitFolder
cd $SubmitFolder

for ((icn=1; icn<=CompNodeNum; icn++))
do
	# Create computational node (cnode) folders
	mkdir cnode$icn
	# Copy femera to the current folder: 1-exe, 2-neper ony
	scp -r ../femera-Backup ./cnode$icn/femera-mini-app-master
	rm -r ./cnode$icn/femera-mini-app-master/neper
	scp -r ../femera-Neper/neper ./cnode$icn/femera-mini-app-master/
	# Copy control file to the current folder
	scp ../AutoRunFemeraExample.sh ./cnode$icn/AutoRunFemera.sh
	# Copy CPU 
	scp ../femera-Backup/cpumodel.sh ./cnode$icn/cpumodel.sh
	# Copy Run.sh
	#sed -i "29s/.*/icn=$icn/" ../Run.sh
	sed -i "33s/.*/icn=$icn/" ../Run.sh

	scp ../Run.sh ./cnode$icn/Run.sh
done
cd ..

# Modify AssembleAMP.sh
sed -i "12s/.*/CompNodeNum=$CompNodeNum/" ./AssembleAMP.sh
sed -i "13s/.*/GrainNum=$GrainNum/" ./AssembleAMP.sh
sed -i "14s/.*/dGrainNum=$dGrainNum/" ./AssembleAMP.sh

$python3 modifyrun.py $GrainNum $dGrainNum $CompNodeNum


