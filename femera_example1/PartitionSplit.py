# This file is designed to apply periodic boundary condition for Phase UCP in Femera
#!/usr/bin/env python
import numpy as np
import os, sys
import math
import csv
import time
import h5py

# Get current work directory
pwd = os.getcwd()

print("========== Split GMSH File According to Partitions ===========")

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum)

# Prefix
prefix = sys.argv[2]
# prefix = int(CompNodeNum)
print('prefix is ',prefix, flush=True)

path = pwd

#=========================================================
#                       Read Elements
#=========================================================
AbaFile=path+'/'+prefix+'_elems.inp'

# Find *Element
tic = time.perf_counter()
HexEls = []
with open(AbaFile) as file:
	data = file.readlines()
	iline = 0
	for dataline in data:
		# print('dataline is: ',dataline)
		if dataline.startswith('*Element, type=C3D8'):
			iline = iline + 1
			while not data[iline].startswith('*'):
				HexEls.append(data[iline].strip())
				iline = iline + 1
		iline = iline + 1
nEls = len(HexEls)
print('nEls is ',nEls)
toc = time.perf_counter()
print("Read the Elements in "+str(toc - tic)+" seconds")

# Convert to int
print("Start to convert the Elements variable format from string to int")
tic = time.perf_counter()
for iEl in range(nEls):
		if iEl % 1000000 == 0:
			print('	'+str(iEl/1000000)+" million elements format converted!") 
		HexEls[iEl] = HexEls[iEl].split(', ')
		HexEls[iEl] = [int(HexEls[iEl][0]),int(HexEls[iEl][1]),int(HexEls[iEl][2]),int(HexEls[iEl][3]),int(HexEls[iEl][4]),int(HexEls[iEl][5]),int(HexEls[iEl][6]),int(HexEls[iEl][7]),int(HexEls[iEl][8])]
toc = time.perf_counter()
HexEls = np.array(HexEls)
print("Convert the HexEls variable format in "+str(toc - tic)+" seconds")

#=========================================================
#                       Read ElSets
#=========================================================
AbaFile=path+'/'+prefix+'_elset.inp'
# Find *Elset, elset=Grain
tic = time.perf_counter()
ElSets = []

with open(AbaFile) as file:
	data = file.readlines()
	iline = 0
	ig = 0
	for dataline in data:
		# print('dataline is: ',dataline)
		if dataline == '*Elset, elset=Grain'+str(ig+1)+'_set\n':
			tmpElSets = []
			ilineOld = iline
			iline = iline + 1
			
			while not data[iline].startswith('*'):
				#print('data[iline] is ',data[iline])
				tmpElSets.append(data[iline].strip())
				iline = iline + 1
			ElSets.append(' '.join(tmpElSets))
			iline = ilineOld
			ig = ig + 1
		iline = iline + 1

nElSets = len(ElSets)
print('nElSets is ',nElSets)
# Convert El IDs to int
for iElSet in range(nElSets):
	ElSets[iElSet] = ElSets[iElSet].replace(',','')
	ElSets[iElSet] = ElSets[iElSet].split(' ')
	for iEl in range(len(ElSets[iElSet])):
		ElSets[iElSet][iEl] = int(ElSets[iElSet][iEl])
for ig in range(GrainNum):
	ElSets[ig] = np.array(ElSets[ig])
#print('ElSets is ',ElSets)
toc = time.perf_counter()
print("Read the ElSets in "+str(toc - tic)+" seconds")

#=========================================================
#                   Element Set IDs
#=========================================================
tic = time.perf_counter()
SetIDs = [0 for x in range(nEls)]
SetIDs = np.array(SetIDs)
for ig in range(GrainNum):
	ElSets[ig] = np.array(ElSets[ig])
	#print('ElSets[ig] is ',ElSets[ig])
	SetIDs[ElSets[ig]-1] = [ig+1]*len(ElSets[ig])
SetIDs = list(SetIDs)
toc = time.perf_counter()
print("Find element set IDs in "+str(toc - tic)+" seconds")
#del globals()['ElSets']

#=========================================================
#             Split Hex and Print Grain Elements
#=========================================================
print("Start to save GrainEls to files")
tic = time.perf_counter()
for ig in range(GrainNum):
	if ig % 1000 == 0:
			print('	'+str(ig/1000)+" thousand grains converted!") 
	GrainEls = HexEls[ElSets[ig]-1]
	# Split Hex element
	TetEls = np.zeros((len(GrainEls),6,10))
	for iEl in range(len(GrainEls)):
		# 1st Tet: 1-2-4-5
		TetEls[iEl,0,:] = [(GrainEls[iEl,0]-1)*6+1, 4, 3, SetIDs[GrainEls[iEl,0]-1], SetIDs[GrainEls[iEl,0]-1], 0, GrainEls[iEl,1], GrainEls[iEl,2], GrainEls[iEl,4], GrainEls[iEl,5] ]
		# 2nd Tet: 2-4-5-8
		TetEls[iEl,1,:] = [(GrainEls[iEl,0]-1)*6+2, 4, 3, SetIDs[GrainEls[iEl,0]-1], SetIDs[GrainEls[iEl,0]-1], 0, GrainEls[iEl,2], GrainEls[iEl,4], GrainEls[iEl,5], GrainEls[iEl,8] ]
		# 3rd Tet: 2-8-5-6
		TetEls[iEl,2,:] = [(GrainEls[iEl,0]-1)*6+3, 4, 3, SetIDs[GrainEls[iEl,0]-1], SetIDs[GrainEls[iEl,0]-1], 0, GrainEls[iEl,2], GrainEls[iEl,8], GrainEls[iEl,5], GrainEls[iEl,6] ]
		# 4th Tet: 2-3-4-8
		TetEls[iEl,3,:] = [(GrainEls[iEl,0]-1)*6+4, 4, 3, SetIDs[GrainEls[iEl,0]-1], SetIDs[GrainEls[iEl,0]-1], 0, GrainEls[iEl,2], GrainEls[iEl,3], GrainEls[iEl,4], GrainEls[iEl,8] ]
		# 5th Tet: 2-3-8-6
		TetEls[iEl,4,:] = [(GrainEls[iEl,0]-1)*6+5, 4, 3, SetIDs[GrainEls[iEl,0]-1], SetIDs[GrainEls[iEl,0]-1], 0, GrainEls[iEl,2], GrainEls[iEl,3], GrainEls[iEl,8], GrainEls[iEl,6] ]
		# 6th Tet: 3-8-6-7
		TetEls[iEl,5,:] = [(GrainEls[iEl,0]-1)*6+6, 4, 3, SetIDs[GrainEls[iEl,0]-1], SetIDs[GrainEls[iEl,0]-1], 0, GrainEls[iEl,3], GrainEls[iEl,8], GrainEls[iEl,6], GrainEls[iEl,7] ]
	h5f = h5py.File(path+'/Partitions/GrainEls/GrainTet'+str(ig+1)+'Els.h5','w')
	h5f.create_dataset('dataset_1',data=TetEls)
	h5f.close()
toc = time.perf_counter()
print("Save GrainEls in "+str(toc - tic)+" seconds")
del globals()['HexEls']

#=========================================================
#                       Read Nodes
#=========================================================
AbaFile=path+'/'+prefix+'_nodes.inp'
# Find *Node
tic = time.perf_counter()
Nodes = []
with open(AbaFile) as file:
	data = file.readlines()
	iline = 0
	for dataline in data:
		# print('dataline is: ',dataline)
		if dataline.startswith('*Node'):
			iline = iline + 1
			while not data[iline].startswith('*'):
				Nodes.append(data[iline].strip())
				iline = iline + 1
		iline = iline + 1
nNodes = len(Nodes)
print('nNodes is ',nNodes)
for iNode in range(nNodes):
	Nodes[iNode] = Nodes[iNode].replace(',','')
toc = time.perf_counter()
print("Read the Nodes in "+str(toc - tic)+" seconds")

# Convert to int
print("Start to convert the Nodes variable format from string to int")
tic = time.perf_counter()
for iNode in range(nNodes):
		if iNode % 1000000 == 0:
			print('	'+str(iNode/1000000)+" million nodes format converted!") 
		Nodes[iNode] = Nodes[iNode].split(' ')
		Nodes[iNode] = [int(Nodes[iNode][0]),float(Nodes[iNode][1]),float(Nodes[iNode][2]),float(Nodes[iNode][3])]
Nodes = np.array(Nodes)
toc = time.perf_counter()
print("Convert the Nodes variable format in "+str(toc - tic)+" seconds")

#=========================================================
#                       Find Grain Nodes
#=========================================================
print("Start to find grain nodes and write new files")
tic = time.perf_counter()
for ig in range(GrainNum):
	if ig % 1000 == 0:
			print('	'+str(ig/1000)+" thousand grains split!") 
	# Read in GrainEls
	h5f = h5py.File(path+'/Partitions/GrainEls/GrainTet'+str(ig+1)+'Els.h5','r')
	GrainEls = h5f['dataset_1'][:]
	h5f.close()
	GrainNodeIDs = GrainEls[:,:,6:10]
#	print('GrainEls is ',GrainEls)
#	print('GrainNodeIDs is ',GrainNodeIDs)
	GrainNodeIDs = GrainNodeIDs.reshape(-1)
	GrainNodeIDs = np.unique(GrainNodeIDs)
	GrainNodeIDs = GrainNodeIDs.astype(int)
#	print('GrainNodeIDs is ',GrainNodeIDs)
	GrainNodes = Nodes[np.array(GrainNodeIDs)-1,:]
	# Write new file
	with open(path+'/Partitions/PartitionSplit_'+str(ig+1)+'.msh',"w") as file:
		file.writelines('$MeshFormat\n2.2 0 8\n$EndMeshFormat\n')
		file.writelines('$Nodes\n')
		file.writelines(str(len(GrainNodes))+'\n')
		np.savetxt(file, GrainNodes, fmt='%i %5f %5f %5f')
		file.writelines('$EndNodes\n')
		file.writelines('$Elements\n')
		file.writelines(str(len(GrainEls)*6)+'\n')
		np.savetxt(file, GrainEls.transpose(0,1,2).reshape(-1,GrainEls.shape[2]), fmt='%5d', delimiter=' ')
		file.writelines('$EndElements\n')
del globals()['Nodes']
toc = time.perf_counter()
print("Extract grain nodes and print new file in "+str(toc - tic)+" seconds")




