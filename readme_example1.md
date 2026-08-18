==================================================
# Overview of the folders
generate_geometry_example1 -> here we generate microstructure geometry using dream3d and extract mesh, neighbor, orientation information
femera_example1 -> here we solve the influence function problems using femera (ROM construction)
rom_simulation_example1 -> here we solve the ROM using abaqus/umat

==================================================
# Generate Geometry

Input files are created using DREAM3D. You can open example1.json with DREAM3D to recover the pipeline used for geometry generation.
If you want to rerun the pipeline, change the directories in appropriate filters to the current working directory.

Before proceeding to ROM construction, generate input files (without submitting job to cluster) :
1. run generate.sh script (on local computer). Input is grain number (change it inside the script) and the prefix of the dream3d files ('example1' in my case)
    this will generate cubic*.ori file with bunge euler angles in radians, 
    extract neighbor information creating Neighbors.txt, GNMap.dat, and NBgrain.dat files, 
    Texture.txti files and copies them to both femera folder.

==================================================
# Femera - solve the IFPs

To construct ROM use following commands:
1. Open example1_nodes.inp file. Scroll down and delete the last node. For some reason, Dream3d adds this node denoted as "999999, 0,0,0" at the end of node list. 
2. Copy femera_example1 folder to cluster
3. Run "qsub PerformanceTest.sh" to submit job to cluster (note that we use Sun Grid Engine(SGE) as the cluster job manager). This file will convert mesh input files from Dream3d (abaqus input file format) into femera inputs as well as assign elastic material properties, orientations, and periodic boundary conditions.
    In this file (PerformanceTest.sh), you should set number of grains, set number of compute nodes, set number of cores per node (for openmp library used by femera), set elastic cubic constants, set voxel resolution used in dream3d, and dimensions of the RVE (length in each direction). 
    relative tolerance for femera is set to 1e-10 but it can be changed in AutoRunFemeraExample.sh file.
    Make sure you set the python path correctly in ./PerformanceTest.sh, AssembleAMP.sh, and femera-Backup/Example1.sh files.
    Output is written to ./log.txt file (you can change this name by modifying SGE command) and ./femera-Neper/Example1Result.txt file.
4. Run "./Submit.sh" which will submit femera jobs to the cluster nodes. Note that number of submissions is equal to number of compute nodes (one for this example).
    modify the file before submitting so that jobs are submitted to appropriate nodes. I could have automated it but some of our nodes are not available most of the time so it is done manually.
    Output is written to ./Jobs/cnode1/log.txt (you can change this name by modifying SGE command) file and ./Jobs/cnode1/femera-mini-app-master/TestResult.txt file.
5. Run "qsub AssembleAMP.sh" to submit job to cluster. This script will gather information solved by femera and assemble it into CoefTens.dat and SparseCoefTens.dat files.
    Output is written to ./amp.txt file (you can change this name by modifying SGE command).

Output will be CoefTens.dat (full coefficient tensors without sparsity) and SparseCoefTens.dat (sparse coefficient tensors) files.

==================================================
# Run Abaqus UMAT simulation

Go to rom_esireport folder
1. Before running ROM, we need more input files. Copy ph1.dat, ph2.dat, GNMap.dat, NBgrain.dat, and Texture.txti from generate_geometry_example1 folder to ./rom_example1/Simulation/ folder. Also, from femera folder copy SparseCoefTens.dat file to the ./rom_example1/Simulation/ folder and rename it as CoefTens.dat.
2. Go to RunSparse.sh file change GrainNum to number of grains in the microstructure.
3. To run the simulation run ./rom_example1/Simulation/caller.sh or ./rom_example1/sub.sh folder (this will submit the tension.inp input file to the cluster)
    The log file of the cluster job will be written to runsparselog and simlog files. Also abaqus output is written to tension.odb. Make sure you set the abaqus path (in RunSparse.sh file), python path (in RunSparse.sh file) and intel compiler path script (in inter.sh file) correctly.

==================================================
# Some additional info
ROM simulation assumes that all grains are assigned phase 2 plastic properties. You can change plastic properties in test.xtali file.
You can visualize the DREAM3D generated microstructure by opening example1.xdmf file using Paraview.

## state variable description
there are 71 state variables per grain and they are stored grain by grain (71 sv for grain 1,71 sv for grain 2, ..., 71 sv for grain N). These state variables are
1. 1-6 -> microscale stress components
2. 7-12 -> microscale elastic strain components
3. 25 -> N/A
4. 26 -> sum of absolute gamma dots over all slip systems
5. 27 -> state variable 26 multiplied by time increment
6. 28-39 -> gamma dot (slip rate) for each slip system (12 slip sytems for FCC)
7. 40-45 -> plastic strain components
8. 46-57 -> resolved shear stress per slip systems
9. 58-63 -> plastic strain rate components
10. 64-67 -> fips
11. 68-71 -> sum of absolute gamma dots over a slip plane