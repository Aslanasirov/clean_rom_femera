# Aslan Nasirov 3/3/22
# EXTRACT NEIGHBOR INFORMATION FROM PREFIX.DREAM3D FILE

import h5py
import numpy as np
import os, sys

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum, flush=True)

neighbor_filename = "Neighbors.txt"
NBgrain_filename = "NBgrain.dat"
GNMap_filename = "GNMap.dat"
grainnum = GrainNum

# file1=open(neighbor_filename,'r')

counter=0
numnn = np.zeros((grainnum))
nns = []
with open(neighbor_filename, 'r') as file:
    # print(counter)
    # print(file.read())
    for line in file:
        # print(counter)
        # print(line.rstrip())
        # curline = line.strip()
        curline = line.split()
        # print(len(curline))
        numnn[counter] = len(curline)
        for grain in range(len(curline)):
            nns.append( curline[grain] )
        counter=counter+1




# f.close()

# oris = oris*180/np.pi
# print(str(oris[0,:]))
numnn = numnn.astype(int)
nnzs = len(nns)
counter = 0
with open(NBgrain_filename, 'w') as file:
    pass
    for i in range(grainnum):
        print(i)
        for j in range(numnn[i]):
            file.write(str(nns[counter])+' ')
            counter=counter+1
        for k in range(grainnum - numnn[i]):
            file.write(str(0)+' ')
        file.write('\n')
    
numnn = numnn.astype(int)
nnzs = len(nns)
counter = 0
# with open(GNMap_filename, 'w') as file:
#     file.write(str(grainnum) + '\n')

#     for i in range(grainnum):
#         print(i)
#         file.write(str(numnn[i])+'\n')
#         for j in range(numnn[i]):
#             file.write(str(nns[counter])+' ')
#             counter=counter+1
#         # for k in range(grainnum - numnn[i]):
#         #     file.write(str(0)+' ')
#         file.write('\n')
#     print(i)
#     # if i == 5:
#     #     sys.exit()

with open(GNMap_filename, 'w') as file:
    file.write(str(grainnum) + '\n')

    for i in range(grainnum):
        print(i)
        file.write(str(numnn[i]+1)+'\n')
        flag=1
        if i+1 < int(nns[counter]):
            file.write(str(i+1)+' ')
            flag=0
        for j in range(numnn[i]):
            
            if i+1 < int(nns[counter]) and flag==1:
                file.write(str(i+1)+' ')
                flag=0

            file.write(str(nns[counter])+' ')
            counter=counter+1
            
            # file.write(str(nns[counter])+' ')
            # counter=counter+1
                
            pass
        if flag==1:
            file.write(str(i+1)+' ')
            flag=0
        # for k in range(grainnum - numnn[i]):
        #     file.write(str(0)+' ')
        file.write('\n')
    print(i)
    # if i == 5:
    #     sys.exit()

print('total num nn -> nnz is ', np.sum(numnn))