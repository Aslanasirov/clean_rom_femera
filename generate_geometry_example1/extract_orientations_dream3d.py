# -*- coding: utf-8 -*-
"""
Created on Sun Oct 24 16:41:54 2021

@author: anasi

# EXTRACT ORIENTATIONS FROM DREAM3D FILES (PREFIX.DREAM3D FILE)
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

# open dream3d file
# file_name = 'example1.dream3d'
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

# 
data = subgroup['CellFeatureData']  # load the data
# data = subgroup['CellData']  # load the data
# print(data)
print('following are the data keys')
for datakey in data.keys():
    print(datakey)
angles = data['EulerAngles']

angles = np.array(angles)
angles = np.delete(angles,0,0)
print('shape of angles is ', angles.shape)
print('angles are ',angles)

np.savetxt('cubic'+str(GrainNum)+'.ori',angles)

f.close()  # close the file
