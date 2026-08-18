# ====== This file is employed to automatically run ElasticUCP simulation with femera ========

# Define rain number and delta grain number
GrainNum=$1
dGrainNum=$2
GrainID1=$3
GrainID2=$4
CNum=$5

cwd=$(pwd)
CPUMODEL=`$cwd/cpumodel.sh`

cd ./femera-mini-app-master

# Make 
make base-omp gmsh2fmr

# Execute
time1="$(date -u +%s.%N)"
#exec 2>> TestResult.txt && ./femerb-$CPUMODEL-gcc -v3 -c ${CNum} -s1 -p "neper/PartitionSplit" -r1e-1 -g ${dGrainNum} -j ${GrainID1} -k ${GrainID2} > TestResult.txt

exec 2>> TestResult.txt && ./femerb-$CPUMODEL-gcc -v3 -c ${CNum} -s1 -p "neper/PartitionSplit" -r1e-10 -g ${dGrainNum} -j ${GrainID1} -k ${GrainID2} > TestResult.txt

# exec 2>> TestResult.txt && ./femerb-i7-4578U-gcc -v3 -c ${CNum} -s1 -p "neper/cubic${GrainNum}s1p1" -r1e-7 > TestResult.txt
time2="$(date -u +%s.%N)"
elapsed="$(bc <<<"$time2-$time1")"
echo "Execution time of femera is "
# runtime=$((time2-time1))
echo $elapsed
#python ReformatAMP.py $GrainNum $dGrainNum
