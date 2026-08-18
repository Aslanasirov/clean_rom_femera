# This file is designed to write AMP tensors in CoefTens.dat file
#!/usr/bin/env python
import numpy as np
import os, sys
import math
import csv
import time
import h5py
#from sklearn.neighbors import NearestNeighbors

# Get current work directory
pwd = os.getcwd()

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum, flush=True)

# delta Grain number
dGrainNum = sys.argv[2]
dGrainNum = int(dGrainNum)
print('dGrainNum is ',dGrainNum, flush=True)

# Number of computational nodes
CompNodeNum = sys.argv[3]
CompNodeNum = int(CompNodeNum)
print('CompNodeNum is ',CompNodeNum, flush=True)

os.chdir(pwd+'/Jobs/cnode'+str(CompNodeNum))

ifppernode=int(GrainNum/CompNodeNum)
# print(ifppernode)

residual = GrainNum - ifppernode * (CompNodeNum-1)
# print(residual)

submitted = ifppernode * (CompNodeNum-1)
# print(submitted)
#from sh import sed

# os.system("sed -i 35s/.*/GrainID1="+str(submitted+1)+"/ Run.sh")
# os.system("sed -i 36s/.*/GrainID2="+str(GrainNum)+"/ Run.sh")

ifps = []
ifppernode=int(GrainNum/CompNodeNum)
print(ifppernode)
residual = GrainNum - ifppernode*CompNodeNum
print(residual)

for cnode in range(CompNodeNum):
    if residual == 0:
        ifps.append(ifppernode)
    else:
        ifps.append(ifppernode+1)
        residual=residual-1
    # print("sum is ",sum(ifps))
    cnode=cnode+1
    print(ifps)

    os.chdir(pwd+'/Jobs/cnode'+str(cnode))
    if cnode == 1:
        os.system("sed -i 35s/.*/GrainID1="+str(1)+"/ Run.sh")
        os.system("sed -i 36s/.*/GrainID2="+str(sum(ifps[0:cnode-2])+ifps[cnode-1])+"/ Run.sh")
        print("node ",str(cnode), " from ", str(1), " to ", str(sum(ifps[0:cnode-2])+ifps[cnode-1]))
    else:
        os.system("sed -i 35s/.*/GrainID1="+str(sum(ifps[0:cnode-1])+1)+"/ Run.sh")
        os.system("sed -i 36s/.*/GrainID2="+str(sum(ifps[0:cnode-1])+ifps[cnode-1])+"/ Run.sh")
        print("node ",str(cnode), " from ", str(sum(ifps[0:cnode-1])+1), " to ", str(sum(ifps[0:cnode-1])+ifps[cnode-1]))
