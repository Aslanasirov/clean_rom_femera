# Aslan Nasirov 3/3/22
# This script creates ph1.dat and ph2.dat files with phase information 
# (they are not extracted from DREAM3D and it is assumed that all grains belong to phase 2)
# Neper orientation file expected as input

import h5py
import numpy as np
import os, sys

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum, flush=True)

phase1_filename = "ph1.dat"
phase2_filename = "ph2.dat"
phase1num = 0
phase2num = GrainNum

# file1=open(phase1_filename,'r')
phase1 = np.zeros((phase1num,1))

# file2=open(phase2_filename,'r')
phase2 = np.zeros((phase2num,1))

for i in range(phase2num):
    phase2[i] = i+1

# WRITE FIRST PHASE FILE
with open(phase1_filename, 'w') as file:
    pass
# WRITE SECOND PHASE FILE
with open(phase2_filename, 'w') as file:
    for i in range(phase2num):
        file.write(str(i+1) + '\n')