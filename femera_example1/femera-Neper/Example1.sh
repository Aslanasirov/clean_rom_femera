N=$1 ; dN=$2 ; C=$3 ; 

L1=${13}; L2=${14}; L3=${15};

c1=$4 ; c2=$5 ; c3=$6
c4=$7 ; c5=$8 ; c6=$9
c7=${10} ; c8=${11} ; c9=${12}

python3=$HOME/mylibs/python3/bin/python3

CPUMODEL=`./cpumodel.sh`

echo "Copy Periodic Tet Mesh"
cwd=$(pwd)
echo pwd: $cwd
scp $cwd/../PBC/* neper/

# neper -T -reg 1 -morpho graingrowth -oricrysym cubic -ori uniform -domain "cube("$L,$L,$L")" -n $N -for tess,ori,tesr -o "neper/cubic"$N
# neper -V "neper/cubic"$N".tess" -datacellcol id -print "neper/cubic"$N
# neper -M "neper/cubic"$N".tess" -elttype hex -dim 3 -order 1 -o "neper/cubic"$N -format msh,inp
# #neper -M "neper/cubic"$N".tess" -elttype hex -cl 0.001 -dim 3 -order 1 -o "neper/cubic"$N -format msh,inp
# neper -V "neper/cubic"$N".msh" -dataelsetcol id -print "neper/cubic"$N"HexMesh"
# neper -V "neper/cubic"$N".msh" -dataelsetcol id -cameracoo 3.462*0.005:5.770*0.005:4.327*0.005 -print "neper/cubic"$N"HexMeshCoo1"
# neper -V "neper/cubic"$N".msh" -dataelsetcol id -cameracoo -3.462*0.005:5.770*0.005:4.327*0.005 -print "neper/cubic"$N"HexMeshCoo2"
# neper -V "neper/cubic"$N".msh" -dataelsetcol id -cameracoo 3.462*0.005:5.770*0.005:-4.327*0.005 -print "neper/cubic"$N"HexMeshCoo3"
# ========== Copy Orientation File ==========
#scp ../femera-Backup/neper/cubic16.ori ./neper/
#scp ../femera-Backup/neper/cubic10.ori ./neper/
scp ../femera-Backup/neper/cubic${N}.ori ./neper/
# ========== Generate Periodic Mesh =========
echo "Generate Periodic Mesh"
#python ./PeriodicMesh.py $N
#scp ../PeriodicCubic9.msh neper/PeriodicCubic9.msh
#N=$(expr $N - 1)
echo "Periodic Mesh Generated"
#scp neper/PeriodicCubic$N.msh neper/cubic${N}s1p1.msh2
# neper -V "neper/cubic"$N".tess,neper/cubic"$N".msh" -dataelsetcol id -print "neper/cubic"$N"Mesh"
# neper -V "neper/cubic"$N".tess,neper/cubic"$N"s1p1.msh2" -dataelsetcol id -print "neper/cubic"$N"Mesh2"
echo "Start gmsh2fmr"
# ========== Apply Boundary Conditions =========
#gmsh neper/PeriodicCubic$N.msh -part_split -part $N -o neper/PartitionSplit -format msh2 -save
#python ./PartitionSplit.py $N
#python ../PartitionPBCs.py $N $dN
#./gmsh2fmr-$CPUMODEL-gcc -x@0.0 -x0 -x@0.0 -y0 -x@0.0 -z0 -y@0.0 -x0 -y@0.0 -y0 -y@0.0 -z0 -z@0.0 -x0 -z@0.0 -y0 -z@0.0 -z0 -x@$L -x0 -x@$L -y0 -x@$L -z0 -y@$L -x0 -y@$L -y0 -y@$L -z0 -z@$L -x0 -z@$L -y0 -z@$L -z0 -M0 -E136.31e9 -N0.37 -G127.40e9 -v3 -p "neper/PartitionSplit" #"neper/cubic"$N"s1p1" #
./gmsh2fmr-$CPUMODEL-gcc -x@0.0 -x0 -x@0.0 -y0 -x@0.0 -z0 -y@0.0 -x0 -y@0.0 -y0 -y@0.0 -z0 -z@0.0 -x0 -z@0.0 -y0 -z@0.0 -z0 -x@$L1 -x0 -x@$L1 -y0 -x@$L1 -z0 -y@$L2 -x0 -y@$L2 -y0 -y@$L2 -z0 -z@$L3 -x0 -z@$L3 -y0 -z@$L3 -z0 -M0 -E1e9 -N0.1 -B "neper/cubic$N.ori" -v3 -p "neper/PartitionSplit"
# ========== Apply Boundary Conditions =========
#echo "Apply boundary conditions"
#python ../PBCs.py $N $dN
#python ../BCs.py $N $dN
#python ../FixBCnElast.py $N $dN $c1 $c2 $c3 $c4 $c5 $c6 $c7 $c8 $c9
echo "Boundary conditions applied"
