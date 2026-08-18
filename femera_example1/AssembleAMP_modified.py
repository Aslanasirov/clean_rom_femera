# ASLAN THIS IS A MODIFIED ASSEMBLE AMP FILE
# NOTE: THIS ONE WAS USED FOR CHALLENGE PROBLEM
# NOTE: ASSUMES IFPS ARE DISTRIBUTED !!!
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

print("========== Reformat AMP Tensors ===========")

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

# =======================ASLAN ==========================
ifppernode=int(GrainNum/CompNodeNum)
residual = GrainNum - ifppernode * (CompNodeNum-1)
submitted = ifppernode * (CompNodeNum-1)
ifps = []
ifppernode=int(GrainNum/CompNodeNum)
residual = GrainNum - ifppernode*CompNodeNum

for cnode in range(CompNodeNum):
    if residual == 0:
        ifps.append(ifppernode)
    else:
        ifps.append(ifppernode+1)
        residual=residual-1
    cnode=cnode+1
print("ipfs are ",ifps)
# =======================ASLAN ==========================


Mijkl = np.zeros((GrainNum,6,6))
print('Mijkl is defined!', flush=True)
PartitionVol = [0 for i in range(GrainNum)]
print('PartitionVol is defined!', flush=True)

# """
# ========================== Extract Mijkl and Pijkl ==============================
# Loop over each computational nodes
tic = time.perf_counter()
nUCP0=99
for icn in range(CompNodeNum):
#for icn in range(CompNodeNum-1,CompNodeNum):
	# ========================== ASLAN ==============================
	# if icn < CompNodeNum-1:
		# nUCP = int(GrainNum / dGrainNum / CompNodeNum)
		# nUCP = ifps[icn]
		# print(ifps[icn])
		# nUCP0 = nUCP

	# else:
	# 	tmp = int(GrainNum / dGrainNum / CompNodeNum)
	# 	nUCP = int(GrainNum / dGrainNum - tmp * icn)
	# 	nUCP = ifps[icn]
	# ========================== ASLAN ==============================
	nUCP = int(GrainNum / dGrainNum / CompNodeNum)
	nUCP = ifps[icn]
	nUCP0 = nUCP

	print('icn is ',icn)
	print('nUCP is ',nUCP)

	# ========================== ASLAN ==============================
	# sum(ifps[0:cnode-1])
	# if icn*nUCP0 % 100 == 0:
	# 	print('	'+str(icn*nUCP0/100)+" hundred IFP extracted!", flush=True)
	if sum(ifps[0:icn]) % 100 == 0:
		print('	'+str(sum(ifps[0:icn])/100)+" hundred IFP extracted!", flush=True)
	# ========================== ASLAN ==============================

	with open(pwd+"/Jobs/cnode"+str(icn+1)+"/femera-mini-app-master/TestResult.txt") as file:
		data = file.readlines()
		
		# Find line number of the starting of AMP for each IFP
		tic1 = time.perf_counter()
		
		IFP_line = np.array([line.count('Print AMP tensors') for line in data])
		IFP_line = np.nonzero(IFP_line)
		IFP_line = np.array([item+1 for item in IFP_line])
		IFP_line = np.array(IFP_line[0])
		IFP_line = np.append(IFP_line,[len(data)+1], axis=0)
		print("IFP line is ", IFP_line)
		toc1 = time.perf_counter()
		print("Find line number for this IFP in "+str(toc1 - tic1)+" seconds", flush=True)
		
		for iu in range(nUCP):
			# Initialize Pijkl for the current IFP
			Pijkl_IFP = np.zeros((GrainNum,6,6))
			# Extract data for certain IFP
			SelectedIFP = data[IFP_line[iu]:IFP_line[iu+1]-1]
			#print("SelectedIFP[0:20] is ",SelectedIFP[0:20])
			
#			if iu % 10 == 0:
#				print('	'+str(iu/10)+"*10 IFP extracted!") 
			
		#	tic1 = time.perf_counter()
			
			for ig1 in range(GrainNum):

			#	tic2 = time.perf_counter()
				for irow in range(6):
					tmp = SelectedIFP[ig1*(1+6+1+6*dGrainNum)+1+irow].split(',')
					for icol in range(6):
						Mijkl[ig1,irow,icol] = float(tmp[icol])

				PartitionVol[ig1] = float(SelectedIFP[ig1*(1+6+1+6*dGrainNum)+1+6])
			#	toc2 = time.perf_counter()
			#	print("	Extract M in "+str(toc2 - tic2)+" seconds", flush=True)

			#	tic2 = time.perf_counter()
				for ig2 in range(dGrainNum):
					# Read in Pijkl(iph1,iph2)
					for irow in range(6):
						tmp = SelectedIFP[ig1*(1+6+1+6*dGrainNum)+1+6+1+ig2*6+irow].split(',')
						for icol in range(6):
							# index of IFP: icn*nUCP*dGrainNum+iu*dGrainNum+ig2
							Pijkl_IFP[ig1,irow,icol] = float(tmp[icol])
			#	toc2 = time.perf_counter()
			#	print("	Extract P in "+str(toc2 - tic2)+" seconds", flush=True)
				
		#	toc1 = time.perf_counter()
		#	print("Extract M and P in "+str(toc1 - tic1)+" seconds", flush=True)
			
			# Save Pijkl_IFP for the current IFP
		#	tic1 = time.perf_counter()
			# np.save(pwd+'/Jobs/Pijkl_IFP_'+str(icn*nUCP0+iu+1)+'.txt',Pijkl_IFP)
			
			np.save(pwd+'/Jobs/Pijkl_IFP_'+str(sum(ifps[0:icn])+iu+1)+'.txt',Pijkl_IFP)
			print("icn is ",icn, " printing ", sum(ifps[0:icn])+iu+1)
		#	h5f = h5py.File(pwd+'/Jobs/Pijkl_IFP_'+str(icn*nUCP0+iu+1)+'.h5','w')
		#	h5f.create_dataset('dataset_1',data=Pijkl_IFP)
		#	h5f.close()
		#	toc1 = time.perf_counter()
		#	print("Save M and P in "+str(toc1 - tic1)+" seconds", flush=True)

toc = time.perf_counter()
print("Extract Mijkl and Pijkl in "+str(toc - tic)+" seconds", flush=True)

# Save Mijkl
h5f = h5py.File(pwd+'/Jobs/Mijkl.h5','w')
h5f.create_dataset('dataset_1',data=Mijkl)
h5f.close()

# Save PartitionVol
h5f = h5py.File(pwd+'/Jobs/PartitionVol.h5','w')
h5f.create_dataset('dataset_1',data=PartitionVol)
h5f.close()

# """

# =================== Compute Grain Volume Fraction ======================
tic = time.perf_counter()

h5f = h5py.File(pwd+'/Jobs/PartitionVol.h5','r')
PartitionVol = h5f['dataset_1'][:]
h5f.close()

toc = time.perf_counter()
print("Read in PartitionVol in "+str(toc - tic)+" seconds", flush=True)

VolFrac = [0 for i in range(GrainNum)]
for ig in range(GrainNum):
	VolFrac[ig] = PartitionVol[ig] / sum(PartitionVol)
	


# ========================== Compute Aijkl ==============================
print("Start to compute Aijkl", flush=True)
tic = time.perf_counter()

Aijkl = np.zeros((GrainNum,6,6))
for ig1 in range(GrainNum):

#	tic1 = time.perf_counter()
	
	# read in Pijkl_IFP
#	tic2 = time.perf_counter()
	
	Pijkl_IFP = np.load(pwd+'/Jobs/Pijkl_IFP_'+str(ig1+1)+'.txt.npy')
	
#	h5f = h5py.File(pwd+'/Jobs/Pijkl_IFP_'+str(ig1+1)+'.h5','r')
#	Pijkl_IFP = h5f['dataset_1'][:]
#	h5f.close()
	
#	toc2 = time.perf_counter()
#	print("	Read Pijkl_IFP in "+str(toc2 - tic2)+" seconds", flush=True)
	Aijkl = np.add(Aijkl,Pijkl_IFP)
#	toc1 = time.perf_counter()
#	print("	Compute "+str(ig1)+"_th grain in "+str(toc1 - tic1)+" seconds", flush=True)
	
for ig in range(GrainNum):
	for irow in range(6):
		for icol in range(6):
			if irow == icol:
				Aijkl[ig][irow][icol] = 1.0 - Aijkl[ig][irow][icol]
			else:
				Aijkl[ig][irow][icol] = 0.0 - Aijkl[ig][irow][icol]

np.save(pwd+'/Jobs/Aijkl.txt',Aijkl)

toc = time.perf_counter()
print("Compute Aijkl in "+str(toc - tic)+" seconds", flush=True)

# ========================== Save AMP Tensors to CoefTens.dat ==============================
#print("VolFrac is ",VolFrac)

# Read in Mijkl
h5f = h5py.File(pwd+'/Jobs/Mijkl.h5','r')
Mijkl = h5f['dataset_1'][:]
h5f.close()

# Read in Aijkl
Aijkl = np.load(pwd+'/Jobs/Aijkl.txt.npy')

# Assemble A and M
tic = time.perf_counter()

AM_Print = [0 for i in range(GrainNum*6*6+GrainNum*6*6)]; print("am print shape is ",len(AM_Print),"\n")
iv = 0
for ig in range(GrainNum):
	for icol in range(6):
		for irow in range(6):
			AM_Print[iv] = Aijkl[ig][irow][icol]
			iv = iv + 1
for ig in range(GrainNum):
	for icol in range(6):
		for irow in range(6):
			AM_Print[iv] = Mijkl[ig,irow,icol]
			iv = iv + 1

toc = time.perf_counter()
print("Assemble A and M in "+str(toc - tic)+" seconds", flush=True)

# Print A, M and P to new file
print("Start to write AMP to new file", flush=True)
with open(pwd+'/CoefTens.dat','w') as file:
	file.write(str(GrainNum)+',\n')
	for ig in range(GrainNum):
		file.write(str("{:.6e}".format(VolFrac[ig]))+',')
	file.write('\n')
	# ~~~~~~~~~~~~~~~~~~ Print A and M ~~~~~~~~~~~~~~~~~~~~~
	AM_LineNum = int(len(AM_Print) / 8)
	# Print first AM_LineNum-1 lines
	for iline in range(AM_LineNum):
		file.write(str("{:.6e}".format(AM_Print[iline*8+0]))+', '+str("{:.6e}".format(AM_Print[iline*8+1]))+', '+str("{:.6e}".format(AM_Print[iline*8+2]))+', '+str("{:.6e}".format(AM_Print[iline*8+3]))+', '+str("{:.6e}".format(AM_Print[iline*8+4]))+', '+str("{:.6e}".format(AM_Print[iline*8+5]))+', '+str("{:.6e}".format(AM_Print[iline*8+6]))+', '+str("{:.6e}".format(AM_Print[iline*8+7]))+',\n')
	# Print the special line of AM_Print if exists (No. of terms < 8)
	NumRestItems = len(AM_Print)-AM_LineNum*8
	for i in range(NumRestItems):
		file.write(str("{:.6e}".format(AM_Print[AM_LineNum*8+i]))+', ')
			
	# ~~~~~~~~~~~~~~~~~~ Print P ~~~~~~~~~~~~~~~~~~~~~
	dGrain = GrainNum
	for igCol in range(0,GrainNum,dGrain):
		print("	Start to assemble "+str(igCol)+" to "+str(igCol+dGrain-1)+" cols", flush=True)
		tic = time.perf_counter()
		
		Pijkl_Print = np.zeros((GrainNum,dGrain,6,6))
		for igRow in range(GrainNum):
		
		#	tic1 = time.perf_counter()
			
			# Read in Pijkl component
			Pijkl_tmp = np.load(pwd+'/Jobs/Pijkl_IFP_'+str(igRow+1)+'.txt.npy')
			Pijkl_Print[igRow,:,:,:] = Pijkl_tmp[igCol:igCol+dGrain,:,:]
			
		#	toc1 = time.perf_counter()
		#	print("	Read "+str(igRow)+" row in "+str(toc1 - tic1)+" seconds", flush=True)
			
		toc = time.perf_counter()
		print("		Read data in "+str(toc - tic)+" seconds", flush=True)
		
		Pijkl_Print = Pijkl_Print.transpose(0,1,3,2)# micro transpose
		Pijkl_Print = Pijkl_Print.reshape(GrainNum,dGrain,-1)
		Pijkl_Print = Pijkl_Print.transpose(1,0,2)# macro transpose
		Pijkl_Print = Pijkl_Print.flatten()
		
		# Print first line. NumRestItems=0 means that this is a new line
		for i in range(8-NumRestItems):
			file.write(str("{:.6e}".format(Pijkl_Print[i]))+', ')
		file.write('\n')
		
		# Print middle lines
		Pijkl_NumLine = int((6*6*dGrain*GrainNum-(8-NumRestItems))/8)
		for iline in range(Pijkl_NumLine):
			file.write(str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+0]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+1]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+2]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+3]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+4]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+5]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+6]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+7]))+',\n')

		# Print the special line of Pijkl_Print if exists (No. of terms < 8)
		NumRestItems = 36*dGrain*GrainNum-(8-NumRestItems+Pijkl_NumLine*8)
		for i in range( NumRestItems ):
			file.write(str("{:.6e}".format(Pijkl_Print[i+8-NumRestItems+Pijkl_NumLine*8]))+', ')
		file.write('\n')
		
		toc = time.perf_counter()
		print("		Print these Cols in "+str(toc - tic)+" seconds", flush=True)

