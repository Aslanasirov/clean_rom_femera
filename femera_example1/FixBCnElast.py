# This file is designed to apply periodic boundary condition for Phase UCP in Femera
#!/usr/bin/env python
import numpy as np
import os, sys
import math
import csv

# Get current work directory
pwd = os.getcwd()

print("========== Apply Fixed Boundary Conditinos ===========")

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum)
# CubeSize
# CubeSize = 0.005

# delta Grain number
dGrainNum = sys.argv[2]
dGrainNum = int(dGrainNum)
print('dGrainNum is ',dGrainNum)

# Elastic coefficients
c1 = sys.argv[3]
c1 = float(c1)
c2 = sys.argv[4]
c2 = float(c2)
c3 = sys.argv[5]
c3 = float(c3)
c4 = sys.argv[6]
c4 = float(c4)
c5 = sys.argv[7]
c5 = float(c5)
c6 = sys.argv[8]
c6 = float(c6)
c7 = sys.argv[9]
c7 = float(c7)
c8 = sys.argv[10]
c8 = float(c8)
c9 = sys.argv[11]
c9 = float(c9)

L1 = sys.argv[12]
L1 = float(L1)
L2 = sys.argv[13]
L2 = float(L2)
L3 = sys.argv[14]
L3 = float(L3)

def isCorner(NodeCoor,L1,L2,L3):
	# Coordinates of the current node
	xCoor = NodeCoor[1]
	yCoor = NodeCoor[2]
	zCoor = NodeCoor[3]
	
	# Check 8 corners
	ck1 = 1 if (xCoor == 0 and yCoor == 0 and zCoor == 0) else 0
	ck2 = 1 if (xCoor == 0 and yCoor == L2 and zCoor == 0) else 0
	ck3 = 1 if (xCoor == L1 and yCoor == 0 and zCoor == 0) else 0
	ck4 = 1 if (xCoor == L1 and yCoor == L2 and zCoor == 0) else 0
	
	ck5 = 1 if (xCoor == 0 and yCoor == 0 and zCoor == L3) else 0
	ck6 = 1 if (xCoor == 0 and yCoor == L2 and zCoor == L3) else 0
	ck7 = 1 if (xCoor == L1 and yCoor == 0 and zCoor == L3) else 0
	ck8 = 1 if (xCoor == L1 and yCoor == L2 and zCoor == L3) else 0
	
	if (ck1+ck2+ck3+ck4+ck5+ck6+ck7+ck8 >= 1) :
		return True
	else:
		return False
#=========================================================
#                       Read
#=========================================================
# In Femera, BCs are defined in each fmr file.
for igrain in range(1,GrainNum+1):
	# Clear old data
	if 'nSurfNode' in globals():
		del globals()['nSurfNode']
	if 'nCornerNode' in globals():
		del globals()['nCornerNode']
	# Initialize elastic data
	ElasticData = ["" for x in range(4)]
	# Find $BC0
	with open(pwd+'/neper/PartitionSplit_'+str(igrain)+'.fmr') as file:
		data = file.readlines()
		iline = 0
		for dataline in data:
			# print('dataline is: ',dataline)
			if dataline == '$BC0\n':
				#print('This is grain ',igrain)
				nDOF = int(data[iline+1])
				# print('nDOF is ',nDOF)
				# BC0DOFLineNum = iline+1
				# data[BC0DOFLineNum] = str(nDOF * dGrainNum * 6) + '\n'
				dataBeforeBC = data[0:iline+2]
				#print('nDOF * dGrainNum * 6 is ',nDOF * dGrainNum * 6)
				# number of surface node
				nSurfNode = int(nDOF / 3)
				SurfNodeIDs = ["" for x in range(nSurfNode)]
				#print('nSurfNode is ',nSurfNode)
				# find surface node IDs
				for iSurf in range(1,nSurfNode+1):
					idline = data[iline+3+3*iSurf-3]
					# print('idline is ',idline)
					idline = idline.split()
					SurfNodeIDs[iSurf-1] = idline[0]
			# elastic properties
			if dataline == '$Elastic\n':
				ElasticData[2] = dataline
				ElasticData[3] = '9 '+str(c1)+' '+str(c2)+' '+str(c3)+' '+str(c4)+' '+str(c5)+' '+str(c6)+' '+str(c7)+' '+str(c8)+' '+str(c9)+'\n'
			if dataline == '$Orientation\n':
				ElasticData[0] = data[iline+0]
				ElasticData[1] = data[iline+1]
			iline = iline + 1

	# find corner nodes from surface nodes
	if 'nSurfNode' in globals():
		CornerNodeIDs = []
		nCornerNode = 0
		for iSurf in range(1,nSurfNode+1):
			NodeID = int(SurfNodeIDs[iSurf-1])
			with open(pwd+'/neper/PartitionSplit_'+str(igrain)+'.fmr') as file:
				data = file.readlines()
				iline = 0
				for dataline in data:
					if dataline == '$VertCoor\n':
						# print('iline is ',iline)
						# print('NodeID is ',NodeID)
						# print('data[iline+1+NodeID] is ',data[iline+1+NodeID])
						NodeCoor = data[iline+1+NodeID]
						NodeCoor = NodeCoor.split(' ')
						for inode in range(1,5):
							NodeCoor[inode-1] = float(NodeCoor[inode-1])
					iline = iline + 1
			if isCorner(NodeCoor, L1, L2, L3):
				CornerNodeIDs.append(NodeID)
				nCornerNode = nCornerNode + 1
		#print("CornerNodeIDs are",CornerNodeIDs)
		#print('nCornerNode is ',nCornerNode)

	# boundary conditions
	if 'nCornerNode' in globals() and nCornerNode != 0:
		dataBCs = ["" for x in range(nCornerNode*dGrainNum*6*3)]
		for iCornerNode in range(1,nCornerNode+1):
			for idof in range(0,dGrainNum*6*3):
				dataBCs[iCornerNode*dGrainNum*6*3-dGrainNum*6*3+idof] = str(CornerNodeIDs[iCornerNode-1]) + ' ' + str(idof) + '\n'
		# print('dataBCs is ',dataBCs)
	
	# write modified file
	if 'nCornerNode' in globals() and nCornerNode != 0:
		dataBeforeBC[-1] = str(nCornerNode*dGrainNum*6*3) + '\n'
		with open(pwd+'/neper/PartitionSplit_'+str(igrain)+'.fmr','w') as file:
			file.writelines(dataBeforeBC)
		with open(pwd+'/neper/PartitionSplit_'+str(igrain)+'.fmr','a') as file:
			file.writelines(dataBCs)
		with open(pwd+'/neper/PartitionSplit_'+str(igrain)+'.fmr','a') as file:
			file.writelines(ElasticData)
		del globals()['nSurfNode']
		del globals()['nCornerNode']
	elif 'nCornerNode' in globals() and nCornerNode == 0:
		del dataBeforeBC[-1]
		del dataBeforeBC[-1]
		with open(pwd+'/neper/PartitionSplit_'+str(igrain)+'.fmr','w') as file:
			file.writelines(dataBeforeBC)
		with open(pwd+'/neper/PartitionSplit_'+str(igrain)+'.fmr','a') as file:
			file.writelines(ElasticData)
		del globals()['nSurfNode']
		del globals()['nCornerNode']
	else:
		#print('Test $Elastic')
		# os.system("sed -i '/Orientation/,+1 d' "+pwd+"/neper/PartitionSplit_"+str(igrain)+".fmr")
		os.system("sed -i '/$Elastic/{n;s/.*/9 "+str(c1)+" "+str(c2)+" "+str(c3)+" "+str(c4)+" "+str(c5)+" "+str(c6)+" "+str(c7)+" "+str(c8)+" "+str(c9)+" "+"/}' "+pwd+"/neper/PartitionSplit_"+str(igrain)+".fmr")






