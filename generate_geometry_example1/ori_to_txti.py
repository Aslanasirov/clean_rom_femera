# Aslan Nasirov 3/3/22
# This script converts orientation units from degrees to radians
# Neper orientation file expected as input

import numpy as np 
import os, sys

# Grain number
GrainNum = sys.argv[1]
GrainNum = int(GrainNum)
print('GrainNum is ',GrainNum, flush=True)

input_filename = "cubic" + str(GrainNum)+".ori"
output_filename = "Texture.txti"


# OPEN INPUT FILE
f=open(input_filename,'r')

oris = np.zeros((GrainNum,3))

count=0
for line in f:
    a = line.split()
    for i in range(3):
        oris[count,i] = a[i]
    count = count +1

print("orientations after transformation \n", oris)
f.close()

# CONVERT TO DEGREES
oris = oris*180/np.pi

# WRITE OUTPUT FILE
with open(output_filename, 'w') as file:
    file.write(str(GrainNum)+'\n')
    file.write(str(0)+' ' + str(0)+'\n')
    print("orientations ", oris)
    np.savetxt(file, oris, delimiter=' ')
    file.write(str(GrainNum)+'\n')