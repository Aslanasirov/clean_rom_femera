# This file is designed to apply periodic boundary condition for Phase UCP in Femera
#!/usr/bin/env python
import numpy as np
import os, sys
import math
import csv
import time
#from sklearn.neighbors import NearestNeighbors

print("========== Apply Periodic Boundary Conditinos in Gmsh ===========")

path = os.getcwd()

#path = '/home/xiaoyu/MyResearch/NASA/ChallengeProblem/Simplified/test/14Grains'

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum)

# CubeSize: (L1, L2, L3)
L1 = sys.argv[2]
L1 = float(L1)
L2 = sys.argv[3]
L2 = float(L2)
L3 = sys.argv[4]
L3 = float(L3)
# Voxel Resolution
R = sys.argv[5]
R = int(R)
print('Voxel Resolution is ',R)

# Prefix
prefix = sys.argv[6]
# prefix = int(CompNodeNum)
print('prefix is ',prefix, flush=True)

# ============================ Define some useful functions ========================
def isCorner(NodeCoor,L1,L2,L3):
	# Coordinates of the current node
	xCoor = NodeCoor[1]
	yCoor = NodeCoor[2]
	zCoor = NodeCoor[3]
	#	print('xCoor is ',xCoor)
	#	print('type of xCoor is ',type(xCoor))
	#	print('xCoor == 0.0 is ',xCoor == 0.0)
	#	exit()
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

def isOnEdge(NodeCoor,L1,L2,L3):
	# Coordinates of the current node
	xCoor = NodeCoor[1]
	yCoor = NodeCoor[2]
	zCoor = NodeCoor[3]
	
	# Check 12 edges
	ck1 = 1 if (xCoor == 0 and yCoor == 0) else 0
	ck2 = 1 if (xCoor == 0 and zCoor == 0) else 0
	ck3 = 1 if (yCoor == 0 and zCoor == 0) else 0
	
	ck4 = 1 if (xCoor == 0 and yCoor == L2) else 0
	ck5 = 1 if (xCoor == 0 and zCoor == L3) else 0
	ck6 = 1 if (yCoor == 0 and zCoor == L3) else 0
	
	ck7 = 1 if (xCoor == L1 and yCoor == 0) else 0
	ck8 = 1 if (xCoor == L1 and zCoor == 0) else 0
	ck9 = 1 if (yCoor == L2 and zCoor == 0) else 0
	
	ck10 = 1 if (xCoor == L1 and yCoor == L2) else 0
	ck11 = 1 if (xCoor == L1 and zCoor == L3) else 0
	ck12 = 1 if (yCoor == L2 and zCoor == L3) else 0
	
	if (ck1+ck2+ck3+ck4+ck5+ck6+ck7+ck8+ck9+ck10+ck11+ck12 >= 1) :
		return True
	else:
		return False

def CountSlave(PeriodicNodes,L1,L2,L3):
	SlaveNodes = []
	for inode in PeriodicNodes:
		xCoor = inode[1]
		yCoor = inode[2]
		zCoor = inode[3]
		
		# Check X,Y,Z surface and Xz,XZ,xZ, Yz,YZ,yZ, Xy, XY,xY edges
		ck1 = 1 if (xCoor == L1) else 0
		ck2 = 1 if (yCoor == L2) else 0
		ck3 = 1 if (zCoor == L3) else 0
		
		ck4 = 1 if (xCoor == L1 and zCoor == 0) else 0
		ck5 = 1 if (xCoor == L1 and zCoor == L3) else 0
		ck6 = 1 if (xCoor == 0 and zCoor == L3) else 0
		
		ck7 = 1 if (yCoor == L2 and zCoor == 0) else 0
		ck8 = 1 if (yCoor == L2 and zCoor == L3) else 0
		ck9 = 1 if (yCoor == 0 and zCoor == L3) else 0
		
		ck10 = 1 if (xCoor == L1 and yCoor == 0) else 0
		ck11 = 1 if (xCoor == L1 and yCoor == L2) else 0
		ck12 = 1 if (xCoor == 0 and yCoor == L2) else 0
		
		if (ck1+ck2+ck3+ck4+ck5+ck6+ck7+ck8+ck9+ck10+ck11+ck12 >= 1)  and not isCorner(inode,L1,L2,L3):
			SlaveNodes.append(inode)
		
	return SlaveNodes

def CountMaster(PeriodicNodes,L1,L2,L3):
	MasterNodes = []
	for inode in PeriodicNodes:
		xCoor = inode[1]
		yCoor = inode[2]
		zCoor = inode[3]
		
		# Check x,y,z surface and xz, yz, xy edges
		ck1 = 1 if (xCoor == 0 and (yCoor != L2) and (zCoor != L3)) else 0
		ck2 = 1 if (yCoor == 0 and (xCoor != L1) and (zCoor != L3)) else 0
		ck3 = 1 if (zCoor == 0 and (xCoor != L1) and (yCoor != L2)) else 0
		
	#	if (ck1+ck2+ck3+ck4+ck5+ck6 >= 1) and not isCorner(inode,L1,L2,L3):
		if (ck1+ck2+ck3 >= 1) and not isCorner(inode,L1,L2,L3):
			MasterNodes.append(inode)
		
	return MasterNodes

def isOnXxSurf(inode,L1,L2,L3):
	if inode[1] == 0 or inode[1] == L1:
		return True
	else:
		return False

def isOnYySurf(inode,L1,L2,L3):
	if inode[2] == 0 or inode[2] == L2:
		return True
	else:
		return False

def isOnZzSurf(inode,L1,L2,L3):
	if inode[3] == 0 or inode[3] == L3:
		return True
	else:
		return False

def isOnXZEdges(inode,L1,L2,L3):
	if (inode[1]==L1 and inode[3]==0) or (inode[1]==0 and inode[3]==L3) or (inode[1]==L1 and inode[3]==L3):
		return True
	else:
		return False

def isOnYZEdges(inode,L1,L2,L3):
	if (inode[2]==L2 and inode[3]==0) or (inode[2]==0 and inode[3]==L3) or (inode[2]==L2 and inode[3]==L3):
		return True
	else:
		return False

def isOnXYEdges(inode,L1,L2,L3):
	if (inode[1]==L1 and inode[2]==0) or (inode[1]==0 and inode[2]==L2) or (inode[1]==L1 and inode[2]==L2):
		return True
	else:
		return False

def PairMasterSlaveEdges(SlaveNodes,MasterNodes):
	MSpairEdges = ["" for x in range(len(SlaveNodes))]
	ind = 0
	for inode in SlaveNodes:
	#	if ind % 1000 == 0:
	#		print('	'+str(ind/1000)+" thousand nodes paired!") 
		if isOnXZEdges(inode,L1,L2,L3):
			SubMaster = [x for x in MasterNodes if (x[2] == inode[2])]
		elif isOnYZEdges(inode,L1,L2,L3):
			SubMaster = [x for x in MasterNodes if (x[1] == inode[1])]
		elif isOnXYEdges(inode,L1,L2,L3):
			SubMaster = [x for x in MasterNodes if (x[3] == inode[3])]
		else:
			print('failed node coordinates are ',inode)
			sys.exit("This node is not edge slave node!")
		MSpairEdges[ind] = [inode[0],SubMaster[0][0]]
		ind = ind + 1
	return MSpairEdges

def PairMS4XxSurfs(SlaveNodes,MasterNodes):
	# Sort the slave and master nodes
	SlaveNodes = np.array(SlaveNodes)
	MasterNodes = np.array(MasterNodes)
	#print('SlaveNodes is ',SlaveNodes)
	SlaveIndex = np.lexsort((SlaveNodes[:,3],SlaveNodes[:,2]))
	SlaveSorted = SlaveNodes[SlaveIndex]
	MasterSorted = MasterNodes[SlaveIndex]
	# Pair
	MSpairXxSurf = ["" for x in range(len(SlaveNodes))]
	ind = 0
	for inode in SlaveSorted:
	#	if ind % 10000 == 0:
	#		print('		'+str(ind/10000)+" ten thousand nodes paired!") 
		SubMaster = MasterSorted[ind,:]
		#print('SubMaster is ',SubMaster)
		MSpairXxSurf[ind] = [inode[0],SubMaster[0]]
		#print('MSpairXxSurf[ind] is ',MSpairXxSurf[ind])
		ind = ind + 1
	return MSpairXxSurf

def PairMS4YySurfs(SlaveNodes,MasterNodes):
	# Sort the slave and master nodes
	SlaveNodes = np.array(SlaveNodes)
	MasterNodes = np.array(MasterNodes)
	SlaveIndex = np.lexsort((SlaveNodes[:,3],SlaveNodes[:,1]))
	SlaveSorted = SlaveNodes[SlaveIndex]
	MasterSorted = MasterNodes[SlaveIndex]
	# Pair
	MSpairYySurf = ["" for x in range(len(SlaveNodes))]
	ind = 0
	for inode in SlaveSorted:
	#	if ind % 10000 == 0:
	#		print('		'+str(ind/10000)+" ten thousand nodes paired!") 
		SubMaster = MasterSorted[ind,:]
		MSpairYySurf[ind] = [inode[0],SubMaster[0]]
		ind = ind + 1
	return MSpairYySurf

def PairMS4ZzSurfs(SlaveNodes,MasterNodes):
	# Sort the slave and master nodes
	SlaveNodes = np.array(SlaveNodes)
	MasterNodes = np.array(MasterNodes)
	SlaveIndex = np.lexsort((SlaveNodes[:,2],SlaveNodes[:,1]))
	SlaveSorted = SlaveNodes[SlaveIndex]
	MasterSorted = MasterNodes[SlaveIndex]
	# Pair
	MSpairZzSurf = ["" for x in range(len(SlaveNodes))]
	ind = 0
	for inode in SlaveSorted:
	#	if ind % 10000 == 0:
	#		print('		'+str(ind/10000)+" ten thousand nodes paired!") 
		SubMaster = MasterSorted[ind,:]
		MSpairZzSurf[ind] = [inode[0],SubMaster[0]]
		ind = ind + 1
	return MSpairZzSurf

def  CountLines(GrainSlaveNodes,VertCoor):
	SlaveLines = []
	LenVert = len(VertCoor)
	if not GrainSlaveNodes:
		return SlaveLines
	#print("VertCoor is ",VertCoor)
	for iline in range(0,LenVert):
		#print("VertCoor[iline] is ",VertCoor[iline])
		#print("GrainSlaveNodes is ",GrainSlaveNodes)
		if NodeIsInNodeSet(VertCoor[iline],GrainSlaveNodes):
			SlaveLines.append(iline)
	for iline in SlaveLines:
		iline = int(iline)
	#print("SlaveLines is ",SlaveLines)
	#print("GrainSlaveNodes is ",GrainSlaveNodes)
	#print("VertCoor is ",VertCoor)
	return SlaveLines

def NodeIsInNodeSet(inode,VertCoor):
	#print("In NodeIsInNodeSet,VertCoor is ",VertCoor)
	if isinstance(VertCoor,list):
		if isinstance(VertCoor[0],list):
			for iNd in VertCoor:
				#print("inode is ",inode)
				#print("iNd is ",iNd)
				if inode[0] == iNd[0]:
					return True
			return False
		else:
			if inode == VertCoor:
				return True
			else:
				return False
	else:
		sys.exit("Error! VertCoor in NodeIsInNodeSet is not a list!")

# ============================ Read in Nodes ========================
AbaFile=path+'/'+prefix+'_nodes.inp'

print('Start to read in nodes')
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
		if iNode % 5000000 == 0:
			print('	'+str(iNode/5000000)+" * 5 million nodes format converted!") 
		#Nodes[iNode] = Nodes[iNode]
		Nodes[iNode] = np.asarray(Nodes[iNode].split(' '),dtype=float)
		#list(map(float,Nodes[iNode]))
		#Nodes[iNode] = [int(Nodes[iNode][0]),float(Nodes[iNode][1]),float(Nodes[iNode][2]),float(Nodes[iNode][3])]
#Nodes = np.array(Nodes)
toc = time.perf_counter()
print("Convert Nodes in "+str(toc - tic)+" seconds")
#print("Nodes are",Nodes)

# ============================ Find Boundary Nodes ========================
# Nodes on boundaries (surface, edge and corner)
print("Start to find nodes on boundaries")
tic = time.perf_counter()
PeriodicNodes = []
EdgeNodes = []
XxSurfNodes = []
YySurfNodes = []
ZzSurfNodes = []
CornerNodes = []

# Sort Nodes
Nodes = np.array(Nodes)
# ===================== CHANGE????? ===================
SortIndex = np.lexsort((Nodes[:,3],Nodes[:,2],Nodes[:,1]))
SortedNodes = Nodes[SortIndex,:]

# Extract corner nodes
L1=L1/R; L2=L2/R; L3=L3/R;
CornerIndex = np.array([0,L3,(L2+1)*(L3+1)-1-L3,(L2+1)*(L3+1)-1,   L1*(L2+1)*(L3+1),L1*(L2+1)*(L3+1)+L3,(L1+1)*(L2+1)*(L3+1)-1-L3,(L1+1)*(L2+1)*(L3+1)-1])
L1=L1*R; L2=L2*R; L3=L3*R;
CornerIndex = CornerIndex.astype(int)
CornerNodes = SortedNodes[CornerIndex]
#print('CornerNodes is ',CornerNodes)
#print('	CornerIndex is ',CornerIndex)
#print("	length of CornerNodes is ",len(CornerNodes))

# Extract edge nodes
L1=L1/R; L2=L2/R; L3=L3/R;
EdgeIndex = np.array(list(range(1,   int(L3),   1))+
			list(range(int(L3)+1,   (int(L2)+1)*(int(L3)+1)-1-int(L3),   int(L3)+1))+
			list(range((int(L3)+1)*(int(L2)+1),   int(L1)*(int(L2)+1)*(int(L3)+1),   (int(L3)+1)*(int(L2)+1)))+
				list(range(int(L3)*2+1,   (int(L2)+1)*(int(L3)+1)-1-1,   int(L3)+1))+
				list(range((int(L2)+1)*(int(L3)+1)-1-int(L3)+1,   (int(L2)+1)*(int(L3)+1)-1,   1))+
				list(range((int(L2)+1)*(int(L3)+1)-1+(int(L3)+1)*(int(L2)+1),   (int(L1)+1)*(int(L2)+1)*(int(L3)+1)-1+1-1,   (int(L3)+1)*(int(L2)+1)))+
			list(range((int(L2)+1)*(int(L3)+1)-1-int(L3)+(int(L3)+1)*(int(L2)+1),   (int(L1)+1)*(int(L2)+1)*(int(L3)+1)-1-int(L3)+1-1,   (int(L3)+1)*(int(L2)+1)))+
			list(range((int(L1)+1)*(int(L2)+1)*(int(L3)+1)-1-int(L3)+1,   (int(L1)+1)*(int(L2)+1)*(int(L3)+1)-1+1-1,   1))+
			list(range(int(L1)*(int(L2)+1)*(int(L3)+1)+int(L3)+1,   (int(L1)+1)*(int(L2)+1)*(int(L3)+1)-1-int(L3)+1-1,   int(L3)+1))+
				list(range(int(L1)*(int(L2)+1)*(int(L3)+1)+1,   int(L1)*(int(L2)+1)*(int(L3)+1)+int(L3)+1-1,   1))+
				list(range(int(L1)*(int(L2)+1)*(int(L3)+1)+int(L3)+int(L3)+1,   (int(L1)+1)*(int(L2)+1)*(int(L3)+1)-1+1-1,   int(L3)+1))+
				list(range(int(L3)+(int(L3)+1)*(int(L2)+1),   int(L1)*(int(L2)+1)*(int(L3)+1)+int(L3)+1-1,   (int(L3)+1)*(int(L2)+1))) )
#print('EdgeIndex is ',EdgeIndex)
L1=L1*R; L2=L2*R; L3=L3*R;
#print('len of EdgeIndex is ',len(EdgeIndex))
EdgeIndex = EdgeIndex.astype(int)
EdgeNodes = SortedNodes[EdgeIndex]
#print('EdgeIndex is ',EdgeIndex)
#print("	length of EdgeNodes is ",len(EdgeNodes))


L1=L1/R; L2=L2/R; L3=L3/R;
# Extract surface nodes
#	Xx surfaces
XxSurfIndex = np.array(list(range(0,(int(L2)+1)*(int(L3)+1)-1,1))+list(range(int(L1)*(int(L2)+1)*(int(L3)+1),(int(L1)+1)*(int(L2)+1)*(int(L3)+1)-1,1)))
XxSurfIndex = XxSurfIndex.astype(int)
XxSurfIndex = np.setdiff1d(XxSurfIndex,EdgeIndex)
XxSurfIndex = np.setdiff1d(XxSurfIndex,CornerIndex)
XxSurfNodes = SortedNodes[XxSurfIndex]
#print('	XxSurfIndex is ',XxSurfIndex)
#print("	length of XxSurfNodes is ",len(XxSurfNodes))

#	Yy surfaces
YySurfIndex = []
for i in range(int(L1)+1):
	YySurfIndex = YySurfIndex + list(range(0+(int(L3)+1)*(int(L2)+1)*i,   int(L3)+(int(L3)+1)*(int(L2)+1)*i,   1))
for i in range(int(L1)+1):
	YySurfIndex = YySurfIndex + list(range((int(L2)+1)*(int(L3)+1)-1-int(L3)+(int(L3)+1)*(int(L2)+1)*i,   (int(L2)+1)*(int(L3)+1)-1+(int(L3)+1)*(int(L2)+1)*i,   1))
YySurfIndex = np.array(YySurfIndex)
YySurfIndex = YySurfIndex.astype(int)
YySurfIndex = np.setdiff1d(YySurfIndex,EdgeIndex)
YySurfIndex = np.setdiff1d(YySurfIndex,CornerIndex)
YySurfNodes = SortedNodes[YySurfIndex]
#print('	YySurfIndex is ',YySurfIndex)
#print("	length of YySurfNodes is ",len(YySurfNodes))

#	Zz surfaces
ZzSurfIndex = []
for i in range(int(L1)+1):
	ZzSurfIndex = ZzSurfIndex + list(range(0 + (int(L3)+1)*(int(L2)+1)*i,   (int(L2)+1)*(int(L3)+1)-1-int(L3) + (int(L3)+1)*(int(L2)+1)*i,   int(L3)+1))
for i in range(int(L1)+1):
	ZzSurfIndex = ZzSurfIndex + list(range(int(L3) + (int(L3)+1)*(int(L2)+1)*i,   (int(L2)+1)*(int(L3)+1)-1 + (int(L3)+1)*(int(L2)+1)*i,   int(L3)+1))
ZzSurfIndex = np.array(ZzSurfIndex)
ZzSurfIndex = ZzSurfIndex.astype(int)
ZzSurfIndex = np.setdiff1d(ZzSurfIndex,EdgeIndex)
ZzSurfIndex = np.setdiff1d(ZzSurfIndex,CornerIndex)
ZzSurfNodes = SortedNodes[ZzSurfIndex]
#print('	ZzSurfIndex is ',ZzSurfIndex)
#print("	length of ZzSurfNodes is ",len(ZzSurfNodes))
L1=L1*R; L2=L2*R; L3=L3*R;

print("	length of CornerNodes is ",len(CornerNodes))
print("	length of EdgeNodes is ",len(EdgeNodes))
print("	length of XxSurfNodes is ",len(XxSurfNodes))
print("	length of YySurfNodes is ",len(YySurfNodes))
print("	length of ZzSurfNodes is ",len(ZzSurfNodes))
toc = time.perf_counter()
print("Found boundary nodes in "+str(toc - tic)+" seconds")

# ============================ Build master-slave connection ============================
print("Start to build master-slave connection")
tic_MSpair = time.perf_counter()
# Build master-slave connection for edge nodes
print("	Start to pair master-slave for edge nodes")
tic = time.perf_counter()
SlaveNodes = CountSlave(EdgeNodes,L1,L2,L3)
MasterNodes = CountMaster(EdgeNodes,L1,L2,L3)
MSpairEdges = PairMasterSlaveEdges(SlaveNodes,MasterNodes)
toc = time.perf_counter()
print("	Edge saster-slave connection built in "+str(toc - tic)+" seconds")

# Build master-slave connection for surface nodes on X
print("	Start to pair master-slave for surface nodes on X")
tic = time.perf_counter()
SlaveNodes = CountSlave(XxSurfNodes,L1,L2,L3)
print("		length of SlaveNodes is ",len(SlaveNodes))
MasterNodes = CountMaster(XxSurfNodes,L1,L2,L3)
print("		length of MasterNodes is ",len(MasterNodes))
MSpairXxSurf = PairMS4XxSurfs(SlaveNodes,MasterNodes)
toc = time.perf_counter()
print("	X Surface master-slave connection built in "+str(toc - tic)+" seconds")

# Build master-slave connection for surface nodes on Y
print("	Start to pair master-slave for surface nodes on Y")
tic = time.perf_counter()
SlaveNodes = CountSlave(YySurfNodes,L1,L2,L3)
print("		length of SlaveNodes is ",len(SlaveNodes))
MasterNodes = CountMaster(YySurfNodes,L1,L2,L3)
print("		length of MasterNodes is ",len(MasterNodes))
MSpairYySurf = PairMS4YySurfs(SlaveNodes,MasterNodes)
toc = time.perf_counter()
print("	Y Surface master-slave connection built in "+str(toc - tic)+" seconds")

# Build master-slave connection for surface nodes on Z
print("	Start to pair master-slave for surface nodes on Z")
tic = time.perf_counter()
SlaveNodes = CountSlave(ZzSurfNodes,L1,L2,L3)
print("		length of SlaveNodes is ",len(SlaveNodes))
MasterNodes = CountMaster(ZzSurfNodes,L1,L2,L3)
print("		length of MasterNodes is ",len(MasterNodes))
MSpairZzSurf = PairMS4ZzSurfs(SlaveNodes,MasterNodes)
toc = time.perf_counter()
print("	Z Surface master-slave connection built in "+str(toc - tic)+" seconds")

MSpair = np.concatenate((MSpairEdges,MSpairXxSurf,MSpairYySurf,MSpairZzSurf),axis=0)
print("		length of MSpair is ",len(MSpair))

toc_MSpair = time.perf_counter()
print("Built master-slave connection in "+str(toc_MSpair - tic_MSpair)+" seconds")

# Extract slave nodes list from MSpair
SlaveInMSpair = MSpair[:,0]
MasterInMSpair = MSpair[:,1]
#print('SlaveInMSpair[0:100] is ',SlaveInMSpair[0:100])
#print('MasterInMSpair[0:100] is ',MasterInMSpair[0:100])
sorter = np.argsort(SlaveInMSpair)

#=========================================================
#                   Write New msh Files
#=========================================================
print("Start to write new msh files")
tic_print = time.perf_counter()
# In Femera, BCs are defined in each msh file.
for igrain in range(1,GrainNum+1):
	#print('	This is grain ',igrain)
	if igrain % 500 == 0:
		print('	'+str(igrain/500)+" *5 hundred grain updated!") 
	
	tic = time.perf_counter()
	# Initialize elastic data
	ElasticData = ["" for x in range(2)]
	# Find $Nodes
	os.system("scp "+path+'/Partitions/PartitionSplit_'+str(igrain)+'.msh '+path+'/PBC/PartitionSplit_'+str(igrain)+'.msh')
	with open(path+'/PBC/PartitionSplit_'+str(igrain)+'.msh') as file:
		data = file.readlines()
		iline = 0
		for dataline in data:
			# print('dataline is: ',dataline)
			if dataline == '$Nodes\n':
				dataBeforeNodes = data[0:iline]
				dataNodesTitle = data[iline:iline+2]
				nNodes = int(data[iline+1])
				dataNodes = data[iline+2:iline+2+nNodes]
				dataNodesEnd = data[iline+2+nNodes]
			if dataline == '$Elements\n':
				dataElTitle = data[iline]
				nEls = int(data[iline+1])
				# Extract element connection table
				dataEls = data[iline+2:iline+2+nEls]
				dataElsEnd = data[iline+2+nEls]
			iline = iline + 1
	toc = time.perf_counter()
#	print("		Check point 1 in "+str(toc - tic)+" seconds")

	# apply periodic boundary conditions
	tic = time.perf_counter()
	# ~~~~~~~~~~ Identify slave nodes ~~~~~~~~~~ 
	GrainNodes = ["" for x in range(nNodes)]
	iv = 0
	for iVert in dataNodes:
		tmp = dataNodes[iv].split()
		GrainNodes[iv] = [int(float(tmp[0])),float(tmp[1]),float(tmp[2]),float(tmp[3])]
		iv = iv + 1
	GrainSlaveNodes = CountSlave(GrainNodes,L1,L2,L3)
	SelectedGrainSlave = GrainSlaveNodes
	if len(SelectedGrainSlave) != 0:
		#print('SelectedGrainSlave is ',SelectedGrainSlave)
		SelectedGrainSlaveIDs = np.array(SelectedGrainSlave)[:,0]
	else:
		continue

	toc = time.perf_counter()
#	print("		Check point 2 in "+str(toc - tic)+" seconds")
	
	# ~~~~~~~~~~ for each slave node, find corresponding master node ~~~~~~~~~~ 
	tic = time.perf_counter()
	if len(SelectedGrainSlave) != 0:
		#print('len of SlaveInMSpair is ',len(SlaveInMSpair))
		#print('len of SelectedGrainSlaveIDs is ',len(SelectedGrainSlaveIDs))
		GrainSlaveIndex = sorter[np.searchsorted(SlaveInMSpair, SelectedGrainSlaveIDs, sorter=sorter)]
		GrainMaster = MasterInMSpair[GrainSlaveIndex]
	toc = time.perf_counter()
#	print("		Check point 3 in "+str(toc - tic)+" seconds")
	
	# ~~~~~~~~~~ replace node ID in dataNodes ~~~~~~~~~~ 
	tic = time.perf_counter()
	#print("dataNodes was ",dataNodes)
	SlaveLines = CountLines(SelectedGrainSlave,GrainNodes)
	for iline in SlaveLines:
		SID = GrainNodes[iline][0]
		MID = GrainMaster[np.where(SelectedGrainSlaveIDs==SID)]
		#print('new MID is ',MID[0])
		#MSpairtmp = [x[1] for x in MSpair if (x[0] == SID)]
		#MID = MSpairtmp[0]
		#print('old MID is ',MID)
		GrainNodes[iline][0] = int(MID[0])
		#print("In PBC modification section, MID is ",MID)
	SpaceString = " "
	for iv in range(len(dataNodes)):
		dataNodes[iv] = str(GrainNodes[iv][0])+' '+str(GrainNodes[iv][1])+' '+str(GrainNodes[iv][2])+' '+str(GrainNodes[iv][3]) + '\n'
	#print("dataNodes is ",dataNodes)
	toc = time.perf_counter()
#	print("		Check point 4 in "+str(toc - tic)+" seconds")
	
	# ~~~~~~~~~~ replace node ID in dataEls ~~~~~~~~~~ 
	tic = time.perf_counter()
	#print("dataEls was ",dataEls)
	GrainEls = ["" for x in range(nEls)]
	ie = 0
	for iEl in dataEls:
		tmp = dataEls[ie].split()
		GrainEls[ie] = [int(tmp[0]),int(tmp[1]),int(tmp[2]),int(tmp[3]),int(tmp[4]),int(tmp[5]),int(tmp[6]),int(tmp[7]),int(tmp[8]),int(tmp[9])]
		ie = ie + 1
	toc = time.perf_counter()
#	print("			Extract GrainEls in "+str(toc - tic)+" seconds")
	# loop over element
	tic1 = time.perf_counter()
	for ie in range(nEls):
		# loop over each node
		for inode in range(4):
			SID = GrainEls[ie][6+inode]
			if SID in SelectedGrainSlaveIDs:
				MID = GrainMaster[np.where(SelectedGrainSlaveIDs==SID)]
				MID = int(MID[0])
			else:
				MID = SID
		#	MSpairtmp = [x[1] for x in MSpair if (x[0] == SID)]
		#	if 0 == len(MSpairtmp):
		#		MID = SID
		#	else:
		#		MID = MSpairtmp[0]
		#		print('Original MID is ',MID)
		#		exit()
			GrainEls[ie][6+inode] = MID
	SpaceString = " "
	toc = time.perf_counter()
#	print("			Replace ID in "+str(toc - tic1)+" seconds")
	tic1 = time.perf_counter()
	for iv in range(len(dataEls)):
		dataEls[iv] = str(GrainEls[iv][0])+' '+str(GrainEls[iv][1])+' '+str(GrainEls[iv][2])+' '+str(GrainEls[iv][3])+' '+str(GrainEls[iv][4])+' '+str(GrainEls[iv][5])+' '+str(GrainEls[iv][6])+' '+str(GrainEls[iv][7])+' '+str(GrainEls[iv][8])+' '+str(GrainEls[iv][9]) + '\n'
	toc = time.perf_counter()
#	print("			Update dataEls in "+str(toc - tic1)+" seconds")
	
	#print("dataEls is ",dataEls)
	toc = time.perf_counter()
#	print("		Check point 5 in "+str(toc - tic)+" seconds")
	
	tic = time.perf_counter()
	# write modified file
	with open(path+'/PBC/PartitionSplit_'+str(igrain)+'.msh','w') as file:
		file.writelines(dataBeforeNodes)
	with open(path+'/PBC/PartitionSplit_'+str(igrain)+'.msh','a') as file:
		file.writelines(dataNodesTitle)
		file.writelines(dataNodes)
		file.writelines(dataNodesEnd)
		file.writelines(dataElTitle)
		file.writelines(str(nEls)+'\n')
		file.writelines(dataEls)
		file.writelines(dataElsEnd)
	toc = time.perf_counter()
#	print("		Check point 6 in "+str(toc - tic)+" seconds")
	#exit()
toc_print = time.perf_counter()
print("Generate new files in "+str(toc_print - tic_print)+" seconds")


