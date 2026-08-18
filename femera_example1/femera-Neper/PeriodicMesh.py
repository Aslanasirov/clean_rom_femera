# This file is designed to apply periodic boundary condition for Phase UCP in Femera
#!/usr/bin/env python
import numpy as np
import os, sys
import math
import csv

# Get current work directory
pwd = os.getcwd()

print("========== Generate Periodic Mesh ===========")

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum)

#=========================================================
#                       Read
#=========================================================
# In gmsh file, each Hex element is divided into 6 Tet elements

# Find $Elements
with open(pwd+'/neper/cubic'+str(GrainNum)+'.msh') as file:
	data = file.readlines()
	iline = 0
	for dataline in data:
		# print('dataline is: ',dataline)
		if dataline == '$Elements\n':
			nEl = int(data[iline+1])
			print('nEl is ',nEl)
			DataBeforeEl = data[0:iline+2]
			DataBeforeEl[-1] = str(nEl*6)+"\n"
			HexElements = ["" for x in range(nEl)]
			# find edge node IDs
			for iEl in range(1,nEl+1):
				HexElements[iEl-1] = data[iline+1+iEl]
			# find edge nodes from surface nodes
		# elastic properties
		if dataline == '$EndElements\n':
			DataAfterEl = data[iline:]
		iline = iline + 1

#=========================================================
#                       Split Hex
#=========================================================
# Format transformation
for iEl in range(1,nEl+1):
	HexElements[iEl-1] = HexElements[iEl-1].split(' ')
	for iterm in range(1,len(HexElements[iEl-1])+1):
		HexElements[iEl-1][iterm-1] = int(HexElements[iEl-1][iterm-1])
#print("HexElements are",HexElements)
# Tet elements
nTet = 6*nEl
TetElements = ["" for x in range(nTet)]
for iEl in range(1,nEl+1):
	# 1st Tet: 1-2-4-5
	TetElements[(iEl-1)*6+0] = [(HexElements[iEl-1][0]-1)*6+1, 4, HexElements[iEl-1][2], HexElements[iEl-1][3], HexElements[iEl-1][4], HexElements[iEl-1][5], HexElements[iEl-1][5+1], HexElements[iEl-1][5+2], HexElements[iEl-1][5+4], HexElements[iEl-1][5+5] ]
	# 2nd Tet: 2-4-5-8
	TetElements[(iEl-1)*6+1] = [(HexElements[iEl-1][0]-1)*6+2, 4, HexElements[iEl-1][2], HexElements[iEl-1][3], HexElements[iEl-1][4], HexElements[iEl-1][5], HexElements[iEl-1][5+2], HexElements[iEl-1][5+4], HexElements[iEl-1][5+5], HexElements[iEl-1][5+8] ]
	# 3rd Tet: 2-8-5-6
	TetElements[(iEl-1)*6+2] = [(HexElements[iEl-1][0]-1)*6+3, 4, HexElements[iEl-1][2], HexElements[iEl-1][3], HexElements[iEl-1][4], HexElements[iEl-1][5], HexElements[iEl-1][5+2], HexElements[iEl-1][5+8], HexElements[iEl-1][5+5], HexElements[iEl-1][5+6] ]
	# 4th Tet: 2-3-4-8
	TetElements[(iEl-1)*6+3] = [(HexElements[iEl-1][0]-1)*6+4, 4, HexElements[iEl-1][2], HexElements[iEl-1][3], HexElements[iEl-1][4], HexElements[iEl-1][5], HexElements[iEl-1][5+2], HexElements[iEl-1][5+3], HexElements[iEl-1][5+4], HexElements[iEl-1][5+8] ]
	# 5th Tet: 2-3-8-6
	TetElements[(iEl-1)*6+4] = [(HexElements[iEl-1][0]-1)*6+5, 4, HexElements[iEl-1][2], HexElements[iEl-1][3], HexElements[iEl-1][4], HexElements[iEl-1][5], HexElements[iEl-1][5+2], HexElements[iEl-1][5+3], HexElements[iEl-1][5+8], HexElements[iEl-1][5+6] ]
	# 6th Tet: 3-8-6-7
	TetElements[(iEl-1)*6+5] = [(HexElements[iEl-1][0]-1)*6+6, 4, HexElements[iEl-1][2], HexElements[iEl-1][3], HexElements[iEl-1][4], HexElements[iEl-1][5], HexElements[iEl-1][5+3], HexElements[iEl-1][5+8], HexElements[iEl-1][5+6], HexElements[iEl-1][5+7] ]

#=========================================================
#                  Write New msh File
#=========================================================
# write modified file
with open(pwd+'/neper/PeriodicCubic'+str(GrainNum)+'.msh','w') as file:
	file.writelines(DataBeforeEl)
with open(pwd+'/neper/PeriodicCubic'+str(GrainNum)+'.msh','a') as file:
	for item in TetElements:
		lineContent = str(item[0])+' '+str(item[1])+' '+str(item[2])+' '+str(item[3])+' '+str(item[4])+' '+str(item[5])+' '+str(item[6])+' '+str(item[7])+' '+str(item[8])+' '+str(item[9])+"\n"
		file.write(lineContent) 
with open(pwd+'/neper/PeriodicCubic'+str(GrainNum)+'.msh','a') as file:
	file.writelines(DataAfterEl)






