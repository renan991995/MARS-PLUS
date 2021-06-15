#!/bin/bash
#PBS -N MARSIL-paper
#PBS -q workq
#PBS -l nodes=cn17:ppn=1
#PBS -o error.out
#PBS -j oe
#PBS -m a
#PBS -M b03504028@ntu.edu.tw

export OMP_NUM_THREADS=1

##################    https://gitlab.dune-project.org/docker/ci/commit/c738e5e1ff5db13509a934288887335131a4e2ce
#export OMPI_MCA_rmaps_base_oversubscribe=1
#export OMPI_MCA_mpi_yield_when_idle=1

# Shut up OpenMPI warning about not being able to use the InfiniBand interface on the Heidelberg nodes, it's not needed anyway
#export OMPI_MCA_btl_base_warn_component_unused=0
##################

cd $PBS_O_WORKDIR

######################################################################

rm  -r  ./CMakeFiles  ./cmake_install.cmake  ./CMakeCache.txt  ./*.o   ./molecule.log   2> /dev/null
rm ./*mds ./*mol ./mds/*enc ./mds/*synonyms  2> /dev/null
#mycmk='/home/akitainu/bin/cmake-install/cmake-2.6.4/bin/cmake'
#mymk='/home/akitainu/bin/make-install/make-3.81/bin/make'
mycmk='/home/tom61212/bin/cmake-install/cmake-3.15.5/bin/cmake'
mymk='/home/tom61212/bin/make-install/make-4.2/bin/make'
mymklist=CMakeLists.txt

infile=6set_memredu.in
logfile=marslog

######################################################################


echo "Job started at :`date '+%Y-%m-%d %H:%M:%S'`"               >  ./time.txt
time1=$SECONDS

##### cmake #####

echo "======= Start calculating ========"
${mycmk} ./${mymklist}
${mymk}
./MARS ./${infile}  > ./${logfile}
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
