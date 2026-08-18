# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 16:41:54 2021

@author: anasi

# EXTRACT NEIGHBOR INFORMATION FROM PREFIX.DREAM3D FILE
"""
import h5py
import numpy as np
import os, sys

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum, flush=True)
# Prefix
prefix = sys.argv[2]
# prefix = int(CompNodeNum)
print('prefix is ',prefix, flush=True)

file_name = prefix+'.dream3d'
f = h5py.File(file_name, 'r+')     # open the file

# List all groups
print('All groups')
for key in f.keys():  # Names of the groups in HDF5 file.
    print(key)
print()

# Get the HDF5 group
group = f['DataContainers']

# Checkout what keys are inside that group.
print('Single group')
for key in group.keys():
    print(key)
print()
# Get the HDF5 subgroup
print('printing groups')
print(group.keys())
# subgroup = group['IN625InitialConditions']
subgroup = group['SyntheticVolumeDataContainer']

print(subgroup)

# Checkout what keys are inside that subgroup.
print('Subgroup')
for key in subgroup.keys():
    print(key)
# print()

# 
# data = subgroup['forhdf']  # load the data
# for key in data.keys():
#     print(key)
# print()
# neighborlist = data['NeighborList']


# 
data = subgroup['CellFeatureData']  # load the data
# data = subgroup['CellData']  # load the data
# print(data)
print('following are the data keys')
for datakey in data.keys():
    print(datakey)
neighborlist = data['NeighborList']

neighborlist = np.array(neighborlist)

print('neighbor list shape is ')
# neighborlist = np.delete(neighborlist,0,0)
print(neighborlist.shape)
# for i in range(neighborlist.shape[0]):
# 	if i>10:
# 		break;
# 	print(neighborlist[i])
print('neighbors are ')
print(neighborlist)

nneighbor = data['NumNeighbors']

nneighbor = np.array(nneighbor)

nneighbor = np.delete(nneighbor,0,0)
print(nneighbor.shape)
count=0; sum=0;
with open("Neighbors.txt",'w') as file:
    for i in range(nneighbor.shape[0]):
        sum = sum + nneighbor[i]
        curlist=[]
        for j in range(int(nneighbor[i])):
            
            file.write(str(neighborlist[count])+' ')
            count=count+1
            
        file.write("\n")
        
pass
print(sum)

f.close()  # close the file
