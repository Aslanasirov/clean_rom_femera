# This file is designed to create new CoefTens for sparse DDEHM model with mock microstructure of the challenge problem
# 	Inputs:
#		CoefTens.dat 	: Old CoefTens for full DDEHM model
#		Neighbors.txt 	: 1st layer grain neighbors
#	Outputs:
#		SparseCoefTens.dat	:	new CoefTens, strain assumption, 1st layer

#!/usr/bin/env python
import numpy as np
import os, sys
import math
import csv
import time
import h5py

print("==========================================", flush=True)
print("           Create New CoefTens            ", flush=True)
print("==========================================", flush=True)

# Get current work directory
pwd = os.getcwd()

# Delete old oututs
if os.path.isfile(pwd+'/SparseCoefTens.dat'):
	os.remove(pwd+'/SparseCoefTens.dat')

#=========================================================
#                 	 Settings
#=========================================================
tic = time.perf_counter()

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum, flush=True)

# Read in grain neighbor info
with open(pwd+'/Neighbors.txt','r') as file:
	GN = file.readlines()
for i in range(0,GrainNum):
	GN[i] = GN[i].strip()
	GN[i] = GN[i].split(' ')
	for j in range(len(GN[i])):
		GN[i][j] = int(GN[i][j])
	GN[i] = np.array(GN[i])
for ig in range(GrainNum):
	GN[ig].sort()

# Compute Grain Volume Fraction
h5f = h5py.File(pwd+'/Jobs/PartitionVol.h5','r')
PartitionVol = h5f['dataset_1'][:]
h5f.close()
vf = [0 for i in range(GrainNum)]
for ig in range(GrainNum):
	vf[ig] = PartitionVol[ig] / sum(PartitionVol)

toc = time.perf_counter()
print("Reading and processing data in "+str(toc - tic)+" seconds", flush=True)

# Read in Mijkl
h5f = h5py.File(pwd+'/Jobs/Mijkl.h5','r')
Mijkl = h5f['dataset_1'][:]
h5f.close()

# Define Aijkl
Aijkl = np.zeros((GrainNum,6,6))

# Computing interaction flags
tic = time.perf_counter()
print('	Computing interaction flags', flush=True)
interactFlags = np.zeros((GrainNum,GrainNum))

for igRow in range(GrainNum):
	for igCol in GN[igRow]:
		interactFlags[igRow,igCol-1] = 1.0
toc = time.perf_counter()
print("	Computing interaction flags in "+str(toc - tic)+" seconds", flush=True)

#print("interactFlags is ", interactFlags, flush=True)

#=========================================================
#                Computing new CoefTens
#=========================================================
print('Computing new CoefTens', flush=True)
tic = time.perf_counter()

# ======================= Update and Print Pijkl =============================
#with open(pwd+'/SparsePijkl.dat','w') as file:
with open(pwd+'/SparseCoefTens.dat','w') as file:

	Pijkl_Diag = np.zeros((GrainNum,6,6))
	
	dGrain = GrainNum
	
	for igRow in range(0,GrainNum,dGrain):
		print("	Start to process "+str(igRow)+" to "+str(igRow+dGrain-1)+" rows", flush=True)
		tic1 = time.perf_counter()
		
		Pijkl_Rows = np.zeros((dGrain,GrainNum,6,6))
		for iRow in range(dGrain):
			# Read in Pijkl component
			Pijkl_tmp = np.load(pwd+'/Jobs/Pijkl_IFP_'+str(igRow+1+iRow)+'.txt.npy')
			Pijkl_Rows[iRow,:,:,:] = Pijkl_tmp
			
		toc1 = time.perf_counter()
		print("		Read data in "+str(toc1 - tic1)+" seconds", flush=True)
		Pijkl_Rows = Pijkl_Rows.transpose(1,0,3,2)# micro and macro transpose
		
		# Update diagonal terms
		tic1 = time.perf_counter()
		for igDiagonal in range(dGrain):
			# Update using terms above
			for ig in range(igRow+igDiagonal):
				if (interactFlags[ig, igRow+igDiagonal] == 0.0) :
					Pijkl_Rows[igRow+igDiagonal,igDiagonal,:,:] = np.add(Pijkl_Rows[igRow+igDiagonal,igDiagonal,:,:], np.multiply(Pijkl_Rows[ig,igDiagonal,:,:], vf[ig]/vf[igRow+igDiagonal]) )
					#Pijkl_Rows[igRow+igDiagonal,igDiagonal,:,:] = [[ Pijkl_Rows[igRow+igDiagonal,igDiagonal,i,j] + vf[ig]/vf[igRow+igDiagonal]*Pijkl_Rows[ig,igDiagonal,i,j] for j in range(6)] for i in range(6)]

			# Update using terms below
			for ig in range(igRow+igDiagonal+1,GrainNum):
				if (interactFlags[ig, igRow+igDiagonal] == 0.0) :
					Pijkl_Rows[igRow+igDiagonal,igDiagonal,:,:] = np.add(Pijkl_Rows[igRow+igDiagonal,igDiagonal,:,:], np.multiply(Pijkl_Rows[ig,igDiagonal,:,:], vf[ig]/vf[igRow+igDiagonal]) )
					#Pijkl_Rows[igRow+igDiagonal,igDiagonal,:,:] = [[ Pijkl_Rows[igRow+igDiagonal,igDiagonal,i,j] + vf[ig]/vf[igRow+igDiagonal]*Pijkl_Rows[ig,igDiagonal,i,j] for j in range(6)] for i in range(6)]
					
		toc1 = time.perf_counter()
		print("         	Update diagonal terms in "+str(toc1 - tic1)+" seconds", flush=True)

		for ig in range(dGrain):
			Pijkl_Diag[igRow+ig,:,:] = Pijkl_Rows[igRow+ig,ig,:,:]
			
		# ~~~~~~~~~~~~~~~~ Update Aijkl ~~~~~~~~~~~~~~~~~~~
		for ig in range(GrainNum):
			for ig1 in range(dGrain):
				if (interactFlags[ig, igRow+ig1] == 1.0) or (ig == igRow+ig1) :
					Aijkl[ig,:,:] = np.add(Aijkl[ig,:,:],Pijkl_Rows[ig,ig1,:,:])

	# Update Aijkl
	for ig in range(GrainNum):
		for irow in range(6):
			for icol in range(6):
				if irow == icol:
					Aijkl[ig,irow,icol] = 1.0 - Aijkl[ig,irow,icol]
				else:
					Aijkl[ig,irow,icol] = 0.0 - Aijkl[ig,irow,icol]

	# ~~~~~~~~~~~~~~ Print A and M ~~~~~~~~~~~~~~~~

	""" Print GrainNum """
	file.write(str(GrainNum)+',\n')
	
	""" Print Volume Fraction """
	for ig in range(GrainNum):
		file.write(str("{:.6e}".format(vf[ig]))+',')
	file.write('\n')
	
	""" Assemble and Print Aijkl and Mijkl """
	# Assemble
	AM_Print = [0 for i in range(GrainNum*6*6+GrainNum*6*6)]
	iv = 0
	for ig in range(GrainNum):
		for icol in range(6): # ASLAN A TRANSPOSED
			for irow in range(6):
				AM_Print[iv] = Aijkl[ig,icol,irow]
				iv = iv + 1
	for ig in range(GrainNum):
		for icol in range(6):
			for irow in range(6):
				AM_Print[iv] = Mijkl[ig,irow,icol]
				iv = iv + 1
	# Print
	AM_LineNum = int(len(AM_Print) / 8)
	for iline in range(AM_LineNum):# Print first AM_LineNum-1 lines
		file.write(str("{:.6e}".format(AM_Print[iline*8+0]))+', '+str("{:.6e}".format(AM_Print[iline*8+1]))+', '+str("{:.6e}".format(AM_Print[iline*8+2]))+', '+str("{:.6e}".format(AM_Print[iline*8+3]))+', '+str("{:.6e}".format(AM_Print[iline*8+4]))+', '+str("{:.6e}".format(AM_Print[iline*8+5]))+', '+str("{:.6e}".format(AM_Print[iline*8+6]))+', '+str("{:.6e}".format(AM_Print[iline*8+7]))+',\n')
	# Print the special line of AM_Print if exists (No. of terms < 8)
	NumRestItems = len(AM_Print)-AM_LineNum*8
	for i in range(NumRestItems):
		file.write(str("{:.6e}".format(AM_Print[AM_LineNum*8+i]))+', ')


	# ~~~~~~~~~~~~~~ Print Pijkl ~~~~~~~~~~~~~~~~

	for igRow in range(0,GrainNum,dGrain):
		print("	Start to print "+str(igRow)+" to "+str(igRow+dGrain-1)+" rows", flush=True)
		# Print diagonal terms and non-zero off-diagonal terms
		tic1 = time.perf_counter()
		OffDiagNum = 0
		for ig in range(dGrain):
			OffDiagNum = OffDiagNum + len(GN[igRow+ig])
		Pijkl_Print = np.zeros((dGrain+OffDiagNum,6,6))
		#print("OffDiagNum is ", OffDiagNum, flush=True)
		
		Pijkl_OffDiag = np.zeros((GrainNum,dGrain,6,6))
		for iRow in range(0,GrainNum,dGrain):
			for ig in range(dGrain):
				# Read in Pijkl component
				Pijkl_tmp = np.load(pwd+'/Jobs/Pijkl_IFP_'+str(iRow+1+ig)+'.txt.npy')
				Pijkl_OffDiag[iRow+ig,0:dGrain,:,:] = Pijkl_tmp[igRow:igRow+dGrain,:,:]
		Pijkl_OffDiag = Pijkl_OffDiag.transpose(0,1,3,2)# micro transpose
		
		iPijkl = 0
		for ig in range(dGrain):
			#print("ig is ", ig, flush=True)
			for ig1 in range(igRow+ig):
				#print("	ig1 is ", ig1, flush=True)
				if (interactFlags[igRow+ig,ig1] == 1.0) :
					#print("		iPijkl is ", iPijkl, flush=True)
					#print("		ig is ", ig, flush=True)
					#print("		ig1 is ", ig1, flush=True)
					Pijkl_Print[iPijkl,:,:] = Pijkl_OffDiag[ig1,ig,:,:]
					iPijkl = iPijkl + 1
			#print("igRow+ig is ",igRow+ig, flush=True)
			#print("iPijkl is ",iPijkl, flush=True)
			#print("Pijkl_Rows[igRow+ig,ig,:,:] is ",Pijkl_Rows[igRow+ig,ig,:,:], flush=True)
			Pijkl_Print[iPijkl,:,:] = Pijkl_Diag[igRow+ig,:,:]
			iPijkl = iPijkl + 1
			#print("iPijkl is ",iPijkl, flush=True)
			for ig1 in range(igRow+ig,GrainNum):
				if (interactFlags[igRow+ig,ig1] == 1.0) :
					Pijkl_Print[iPijkl,:,:] = Pijkl_OffDiag[ig1,ig,:,:]
					iPijkl = iPijkl + 1
					#print("		iPijkl is ",iPijkl, flush=True)
			#print("---- iPijkl is ",iPijkl, flush=True)
		toc1 = time.perf_counter()
		print("		Assemble sparse components in "+str(toc1 - tic1)+" seconds", flush=True)
		#print("		Pijkl_Print[0,:,:] is ", Pijkl_Print[0,:,:], flush=True)

		# --------------- Start to print --------------
		Pijkl_Print = Pijkl_Print.flatten()
		# Print first line. NumRestItems=0 means that this is a new line
		for i in range(8-NumRestItems):
			file.write(str("{:.6e}".format(Pijkl_Print[i]))+', ')
		file.write('\n')
		
		# Print middle lines
		Pijkl_NumLine = int((6*6*(dGrain+OffDiagNum)-(8-NumRestItems))/8)
		for iline in range(Pijkl_NumLine):
			file.write(str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+0]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+1]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+2]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+3]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+4]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+5]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+6]))+', '+str("{:.6e}".format(Pijkl_Print[iline*8+8-NumRestItems+7]))+',\n')

		# Print the special line of Pijkl_Print if exists (No. of terms < 8)
		NumRestItems = 36*(dGrain+OffDiagNum)-(8-NumRestItems+Pijkl_NumLine*8)
		for i in range( NumRestItems ):
			file.write(str("{:.6e}".format(Pijkl_Print[i+8-NumRestItems+Pijkl_NumLine*8]))+', ')
		
		
		toc = time.perf_counter()
		print("		Print these Cols in "+str(toc - tic)+" seconds", flush=True)
		
	file.write('\n')

	#	# Print to file
	#	tic1 = time.perf_counter()
	#	for iline in range(dGrain+OffDiagNum):
	#		for irow in range(6):
	#			for icol in range(6):
	#				file.write(str("{:.6e}".format(Pijkl_Print[iline,irow,icol]))+', ')
	#			file.write('\n')
	#	toc1 = time.perf_counter()
	#	print("		Write to file in "+str(toc1 - tic1)+" seconds", flush=True)

"""

with open(pwd+'/SparseCoefTens.dat','w') as file:

	# Print GrainNum
	file.write(str(GrainNum)+',\n')
	
	# Print Volume Fraction
	for ig in range(GrainNum):
		file.write(str("{:.6e}".format(vf[ig]))+',')
	file.write('\n')
	
	# Assemble and Print Aijkl and Mijkl
	# Assemble
	AM_Print = [0 for i in range(GrainNum*6*6+GrainNum*6*6)]
	iv = 0
	for ig in range(GrainNum):
		for icol in range(6):
			for irow in range(6):
				AM_Print[iv] = Aijkl[ig,irow,icol]
				iv = iv + 1
	for ig in range(GrainNum):
		for icol in range(6):
			for irow in range(6):
				AM_Print[iv] = Mijkl[ig,irow,icol]
				iv = iv + 1
	# Print
	AM_LineNum = int(len(AM_Print) / 8)
	for iline in range(AM_LineNum):# Print first AM_LineNum-1 lines
		file.write(str("{:.6e}".format(AM_Print[iline*8+0]))+', '+str("{:.6e}".format(AM_Print[iline*8+1]))+', '+str("{:.6e}".format(AM_Print[iline*8+2]))+', '+str("{:.6e}".format(AM_Print[iline*8+3]))+', '+str("{:.6e}".format(AM_Print[iline*8+4]))+', '+str("{:.6e}".format(AM_Print[iline*8+5]))+', '+str("{:.6e}".format(AM_Print[iline*8+6]))+', '+str("{:.6e}".format(AM_Print[iline*8+7]))+',\n')
	# Print the special line of AM_Print if exists (No. of terms < 8)
	NumRestItems = len(AM_Print)-AM_LineNum*8
	for i in range(NumRestItems):
		file.write(str("{:.6e}".format(AM_Print[AM_LineNum*8+i]))+', ')

# Merge two files
with open(pwd+'/SparseCoefTens.dat','a+') as file1:
	with open(pwd+'/SparsePijkl.dat','r') as file2:
		file1.write(file2.read())
"""

