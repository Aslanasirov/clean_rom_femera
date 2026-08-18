# This file is designed to apply periodic boundary condition for Phase UCP in Femera
#!/usr/bin/env python
import numpy as np
import os, sys
import math
import csv

# Get current work directory
pwd = os.getcwd()

print("========== Split GMSH File According to Partitions ===========")

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum)

#=========================================================
#                       Read
#=========================================================
# In gmsh file, each Hex element is divided into 6 Tet elements

# Find $Elements
with open(pwd+'/neper/PeriodicCubic'+str(GrainNum)+'.msh') as file:
	data = file.readlines()
	iline = 0
	for dataline in data:
		# Extract Nodes
		if dataline == '$Nodes\n':
			dataBeforeNodes = data[0:iline]
			dataNodesTitle = data[iline]
			nNodes = int(data[iline+1])
			dataNodes = ["" for x in range(nNodes)]
			# find edge node IDs
			for inode in range(nNodes):
				dataNodes[inode] = data[iline+2+inode]
			dataNodesEnd = data[iline+2+nNodes]
		# Extract Elements
		if dataline == '$Elements\n':
			dataElsTitle = data[iline]
			nEls = int(data[iline+1])
			dataElements = ["" for x in range(nEls)]
			# find edge node IDs
			for iEl in range(nEls):
				dataElements[iEl] = data[iline+2+iEl]
			dataElsEnd = data[iline+2+nEls]
		if dataline == '$PhysicalNames\n':
			dataAfterEl = data[iline:]
			nSets = int(data[iline+1])
		iline = iline + 1

#=========================================================
#                       Split Partitions
#=========================================================
# Format transformation
Elements = ["" for x in range(nEls)]
for iEl in range(nEls):
	Elements[iEl] = dataElements[iEl].split(' ')
	for iterm in range(len(Elements[iEl])):
		Elements[iEl][iterm] = int(Elements[iEl][iterm])
Nodes = ["" for x in range(nNodes)]
for iNode in range(nNodes):
	Nodes[iNode] = dataNodes[iNode].split(' ')
	for iterm in range(len(Nodes[iNode])):
		Nodes[iNode][iterm] = float(Nodes[iNode][iterm])
#print("Elements is ", Elements)
#print("Nodes is ", Nodes)

#=========================================================
#                  Write New msh Files
#=========================================================
for igrain in range(1,GrainNum+1):
	# ========== Select elements in the current grain ===========
	GrainElements = []
	GrainElArray = []
	for iEl in range(nEls):
		if Elements[iEl][3] == igrain:
			GrainElements.append(dataElements[iEl])
			GrainElArray.append(Elements[iEl])
	nGrainEls = len(GrainElements)
	print("GrainElements is ",GrainElements)
	print("GrainElArray is ",GrainElArray)
	# ============ Select nodes in the current grain ==============
	# Collect node IDs that belong to the current grain
	GrainNodesArray = []
	for iEl in GrainElArray:
		for inode in iEl[6:10]:
			GrainNodesArray.append(inode)
	GrainNodesArray = np.unique(GrainNodesArray)
	print("GrainNodesArray is ",GrainNodesArray)
	# Collect nodal coordinates
	GrainNodes = []
	for inode in Nodes:
		if inode[0] in GrainNodesArray:
			GrainNodes.append(str(int(inode[0]))+" "+" ".join(repr(e) for e in inode[1:])+'\n')
	print("GrainNodes is ", GrainNodes)
	nGrainNodes = len(GrainNodes)
	# write modified file
	#NewFileName = 'PartitionSplit'+str(GrainNum)+'_'
	NewFileName = 'PartitionSplit_'
	with open(pwd+'/neper/'+NewFileName+str(igrain)+'.msh','w') as file:
		file.writelines(dataBeforeNodes)
	with open(pwd+'/neper/'+NewFileName+str(igrain)+'.msh','a') as file:
		file.writelines(dataNodesTitle)
		file.writelines(str(nGrainNodes)+'\n')
		file.writelines(GrainNodes)
		file.writelines(dataNodesEnd)
		file.writelines(dataElsTitle)
		file.writelines(str(nGrainEls)+'\n')
		file.writelines(GrainElements)
		file.writelines(dataElsEnd)
		#file.writelines(dataAfterEl)
	#os.system("scp "+pwd+'/neper/PartitionSplit_'+str(igrain)+'.msh '+pwd+'/neper/PartitionSplit_NoPBC_'+str(igrain)+'.msh')






