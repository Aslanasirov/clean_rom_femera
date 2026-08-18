#!/bin/bash
# TRANSFER FILES BETWEEN COMPUTER AND CLUSTER

# scp -r /home/xiaoyu/aslan/femera_scaling/* anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/
# scp -r /home/xiaoyu/aslan/femera_scaling/cubic*.ori anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/
# scp -r /home/xiaoyu/aslan/femera_scaling/extract_orientations_dream3d.py anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/

# scp -r /home/xiaoyu/aslan/femera_scaling/mesh_res_10* anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/
# scp -r /home/xiaoyu/aslan/femera_scaling/cubic7961.ori anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/
# scp -r /home/xiaoyu/aslan/femera_scaling/*reference/FixBCnElast.py anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/femer*check/


# scp -r anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/scaling_from_cluster.xdmf /home/xiaoyu/aslan/femera_scaling/femera_res1_paraview/scaling_from_cluster.xdmf
# scp -r anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/IN625InitialConditions.dream3d /home/xiaoyu/aslan/femera_scaling/femera_res1_paraview/
# scp -r anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/scaling_from_cluster.dream3d /home/xiaoyu/aslan/femera_scaling/femera_res1_paraview/

# scp -r /home/xiaoyu/aslan/femera_scaling/*fast/FixBCnElast.py anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/femera_scaling/femera_res1_fixbcspeedcheck/
# scp -r anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/clean_rom_femera/challenge_rom_femera/femera/tarball.tar.gz /home/xiaoyu/aslan/femera_scaling/
# scp -r anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/clean_rom_femera/tarball_parallel.tar.gz /home/xiaoyu/aslan/femera_scaling/

# scp -r ./femera_res10_esireport/* anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/femera_res10_esireport/ 
# scp -r ./femera_res10_esireport/Neighbors.txt anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/femera_res10_esireport/ 
# scp -r ./femera_res10_esireport/Ass*.py anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/femera_res10_esireport/ 

scp -r ./ph1.dat anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/rom_example1/Simulation/
scp -r ./ph2.dat anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/rom_example1/Simulation/
scp -r ./NBgrain.dat anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/rom_example1/Simulation/
scp -r ./GNMap.dat anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/rom_example1/Simulation/
scp -r ./Texture.txti anasirov@cee-mcml-head.vuse.vanderbilt.edu:/home/anasirov/nasa/esireport/rom_example1/Simulation/