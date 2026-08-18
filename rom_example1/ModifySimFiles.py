# This file is designed to modify input files of Abaqus simulation
# IT MODIFIES TEST.XTALI-> NUMBER OF PARTS , TENSION.INP -> NUMBER OF STATE VARIABLES
# IT IS ASSUMED THERE ARE 71 STATE VARIABLES PER PART
# ALSO MODIFED MODULES.F90 WITH NUMBER OF PARTS PER PHASE

#!/usr/bin/env python



import numpy as np
import os, sys
import math
import csv

# Get current work directory
pwd = os.getcwd()

print("========== Modify Simulation Files ===========")

#=========================================================
#                       Reading
#=========================================================
print("~~~ Reading ~~~")

# Read in particle number and cell number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum)

particlenum = GrainNum
NumGrain = GrainNum


# Read in phase info
phase1num = 0
phase2num = NumGrain

#=========================================================
#                       Modules.f90
#=========================================================
print("~~~ Modules.f90 ~~~")
# Open Modules.f90 file
with open(pwd+'/Simulation/Modules.f90') as file:
	data = file.readlines()
	iline = 0
	for dataline in data:
		#print(dataline[0:48])
		# Replace with new parameters
		if dataline[0:48]=='    integer(ikind), parameter   :: NPh=2, NPart=':
			#print('dataline was',dataline)
			dataline=dataline[0:48]+str(NumGrain)+', PhSlip(NPH)=(/12,12/), ParPerPh(NPh)=(/'+str(phase1num)+','+str(phase2num)+'/), maxnumslip=12'+'\n'

			data[iline] = dataline
			#print('dataline is',dataline)
		iline = iline + 1

# write Modules.f90 file
with open(pwd+'/Simulation/Modules.f90','w') as file:
	file.writelines(data)

#=========================================================
#                   test.xtali
#=========================================================
print("~~~ test.xtali ~~~")
# Open test.xtali
with open(pwd+'/Simulation/test.xtali') as file:
	data = file.readlines()
	data[0] = '2    '+str(NumGrain)+'                               / NPH, NPart/'+'\n'
	
# write test.xtali
with open(pwd+'/Simulation/test.xtali','w') as file:
	file.writelines(data)

#=========================================================
#                	   tension.inp
#=========================================================
print("~~~ tension.inp ~~~")
# Open tension.inp
with open(pwd+'/Simulation/tension.inp') as file:
	data = file.readlines()
	iline = 0
	for dataline in data:
		# Replace with new parameters
		if dataline[0:7]=='*Depvar':
			#print 'dataline was',dataline
			#SDVnum = phase2num*459+phase1num*297+NumGrain*2+3+7
			SDVnum = phase2num*71+phase1num*71
			#SDVnum = phase2num*63+phase1num*63
			data[iline+1] = str(SDVnum) +',\n'
			#print 'dataline is',dataline

		# if dataline[0:31] == '*Element Output, directions=YES':
			# data[iline+1] = 'SDV' + str(SDVnum) + ',LE,E,S\n'
		iline = iline + 1

# write tension.inp
with open(pwd+'/Simulation/tension.inp','w') as file:
	file.writelines(data)
#=========================================================
#                Find 2nd Layer Neighbors
#=========================================================
