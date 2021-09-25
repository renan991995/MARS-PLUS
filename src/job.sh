#!/bin/bash
#PBS -N MARS-PLUS
#PBS -q workq
#PBS -l nodes=cn1:ppn=1
#PBS -o error.out
#PBS -j oe
#PBS -m a
#PBS -M b03504028@ntu.edu.tw

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

######################################################################

rm  -r  ./CMakeFiles  ./cmake_install.cmake  ./CMakeCache.txt  ./*.o   ./molecule.log   2> /dev/null
rm ./*mds ./*mol ./mds/*enc ./mds/*synonyms  2> /dev/null
mycmk="/home/tom61212/bin/cmake-install/cmake-3.15.5/bin/cmake"
mymk="/home/tom61212/bin/make-install/make-4.2/bin/make"
mymklist="./CMakeLists.txt"

infile="../INPUTS/control.in"
logfile="./marslog"

######################################################################


echo "Job started at :`date '+%Y-%m-%d %H:%M:%S'`"               >  ./time.txt
time1=$SECONDS

##### cmake #####

echo "======= Start calculating ========"
${mycmk} ./${mymklist}
${mymk}
./MARS ${infile}  > ${logfile}
#mpirun -np 1 -mca btl ^openib ./MARS  ./$infile  >  ./$logfile
#valgrind --track-origins=yes --leak-check=full --trace-children=yes --gen-suppressions=all --log-file=memcheck.log  -v ./MARS  ./$infile  >  ./$logfile
#valgrind --track-origins=yes --leak-check=full -v ./MARS ./$infile  > ./$logfile
echo "======= End calculating ========"



time2=$SECONDS
echo "Job ended at   :`date '+%Y-%m-%d %H:%M:%S'`"              >> ./time.txt
#for calculating elapsed time
s=`echo $((time2-time1))`
sec=`echo $((s%60))`
min1=`echo $((s/60))`
min=`echo $((min1%60))`
hour1=`echo $((min1/60))`
hour=`echo $((hour1%24))`
day=`echo $((hour1/24))`

echo "Running this job took you $day Days $hour Hours $min Mins $sec Secs" >> ./time.txt
