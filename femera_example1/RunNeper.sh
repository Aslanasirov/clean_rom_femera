# ====== This file is employed to automatically run ElasticUCP simulation with femera ========

# Define rain number and delta grain number
GrainNum=$1
dGrainNum=$2
CNum=$3

c1=$4
c2=$5
c3=$6
c4=$7
c5=$8
c6=$9
c7=${10}
c8=${11}
c9=${12}

L1=${13}
L2=${14}
L3=${15}

# Make 
make base-omp gmsh2fmr

# Generate geometry mesh
time1="$(date -u +%s.%N)"
./Example1.sh $GrainNum $dGrainNum $CNum $c1 $c2 $c3 $c4 $c5 $c6 $c7 $c8 $c9 $L1 $L2 $L3>> ./Example1Result.txt
# Execute
time2="$(date -u +%s.%N)"
elapsed="$(bc <<<"$time2-$time1")"
echo "Execution time of Example1.sh is "
# runtime=$((time2-time1))
echo $elapsed

