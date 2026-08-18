# ====== This file is employed to generate input files from DREAM3D ========

# get the current working directory
cwd=$(pwd)
echo pwd: $cwd

# Define parameters
GrainNum=13 # INPUT : NUMBER OF GRAINS
Prefix='example1' # INPUT : DREAM3D PREFIX

# extract orientations
echo ===== extract orientaion =====
python3 extract_orientations_dream3d.py $GrainNum $Prefix
# convert ori to txti
echo ===== convert orientaion to txti =====
python3 ori_to_txti.py $GrainNum
# extract neighbors
echo ===== generate neighbor list =====
python3 extract_neighbors_dream3d.py $GrainNum $Prefix
# generate phase information (here I assume that all grains belong to phase 2)
echo ===== make phases =====
python3 make_phases.py $GrainNum
# generate additional neighbors
echo ===== make additional neighbors =====
python3 make_neighbors.py $GrainNum

echo ===== copy generated files to femera folder=====
scp cubic$GrainNum.ori ../femera_example1/femera-Backup/neper/
scp $Prefix*.inp ../femera_example1/
scp Neighbors.txt ../femera_example1/