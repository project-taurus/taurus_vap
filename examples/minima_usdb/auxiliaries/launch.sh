#!/bin/bash

# For debugging purposes
#set -x

#################
#  Environment  #
#################

inter="usdb.sho"  

omp="no"
nthreads=4

mpi="no"
nprocess=4

heredir=$(pwd)
exedir=$heredir/../../exe
wrkdir=$heredir/wrk
intdir=$heredir/../../extras/hamiltonians

code=taurus_vap.exe

input=input.txt

#################
#  Calculation  #
#################

if [ ! -d $wrkdir ]; then mkdir $wrkdir; fi 

cd $wrkdir

cp $exedir/$code .
cp $intdir/$inter .
cp $heredir/$input .        

# For OpenMP calculations
if [ $omp = "yes" ] && [ $nthreads -ne 0 ]; then
  export OMP_NUM_THREADS=$nthreads
  # only for large model spaces
  #export OMP_STACKSIZE=100M
fi

# Runs the code (MPI or not)
if [ $mpi = "yes" ]; then
  mpirun -np $nprocess ./$code < $input
else
  ./$code < $input
fi 

# Some cleaning
rm -f $code $inter 

