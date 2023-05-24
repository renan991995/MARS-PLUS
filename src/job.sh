#!/bin/bash
#PBS -N MARS-PLUS
#PBS -q workq
#PBS -l nodes=1:ppn=1
#PBS -o MARS-PLUS.out
#PBS -j oe
#PBS -m a

source ~/.bashrc
conda activate MARS+

export OMP_NUM_THREADS=1

######################################################################

cd $PBS_O_WORKDIR
rm  -r  ./CMakeFiles  ./cmake_install.cmake  ./CMakeCache.txt  ./*.o   ./molecule.log   2> /dev/null
cmake CMakeLists.txt 
make -j4

cd $PBS_O_WORKDIR/../LOG_FILES
mkdir ./mds 2> /dev/null
rm ./mds/* 2> /dev/null
#rm ./*mds ./*mol ./mds/*enc ./mds/*synonyms  2> /dev/null

cd $PBS_O_WORKDIR

infile=../INPUTS/control.in
logfile=../LOG_FILES/marslog

######################################################################


echo "Job started at :`date '+%Y-%m-%d %H:%M:%S'`"               >  ./time.txt
time1=$SECONDS

##### cmake #####

echo "======= Start calculating ========"

./MARS-PLUS ${infile}  > ${logfile}
echo "======= End calculating ========"

conda deactivate

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
