#!/bin/bash
#PBS -N SMI_Enumerator
#PBS -q work-8cpu
#PBS -l nodes=cn17:ppn=1
#PBS -o error.out
#PBS -j oe
#PBS -m a

export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

######################################################################

rm  -r  ./CMakeFiles  ./cmake_install.cmake  ./CMakeCache.txt  2>  /dev/null
mycmk='/home/tom61212/bin/cmake-install/cmake-3.15.5/bin/cmake'
mymk='/home/tom61212/bin/make-install/make-4.2/bin/make'
mymklist="./CMakeLists.txt"

infile="./control.in"
logfile="./Enumlog"

######################################################################


echo "Job started at :`date '+%Y-%m-%d %H:%M:%S'`"               >  ./time.txt
time1=$SECONDS




echo "======= Start calculating ========"
${mycmk} ./${mymklist}
${mymk}  
./SMI_Enumerator ${infile} > ${logfile}
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
