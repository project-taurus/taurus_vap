#!/bin/bash

#################
#  Description  #
#################
 
# This is an example of script to compile TAURUS_vap. The code requires the
# BLAS/LAPACK libraries. When using the intel compiler "ifort", we recommend
# to use of their specific Math Kernel Library (MKL).

# The scripts takes two arguments:
#  FC = $1 (Fortran compiler)
#     = "gfortran", "ifort", "mpif90" or "mpiifort"
#
#  TH = $2 (Activate OpenMP threading)
#     = "omp" or "none"

# This script is only given as an example and we do not guarantee that it will
# work on your system. In particular, check the version of your compiler and
# the directories of the libraries.

#################
#  Directories  #
#################

heredir=$(pwd)
srcdir=$heredir/src
wrkdir=$heredir/wrk
exedir=$heredir/exe

#############################################
#  Fortran compiler, options and libraries  #
#############################################

FC=$1

# By default, use gfortran
if [ -z $FC ]; then
 FC="gfortran"
fi   

if [ $FC = "ifort" ] || [ $FC = "mpiifort" ]; then
  LIB=""
  OPT="-O3 -mkl" 
elif [ $FC = "gfortran" ] || [ $FC = "mpif90" ]; then
  LIB="-L/usr/lib -llapack -lblas"
  #OPT="-Wall -Wno-maybe-uninitialized"  # To check the warnings
  OPT="-O3" 
else
  echo "Wrong compiler ('gfortran', 'ifort', 'mpif90' or 'mpiifort'). Exiting."
  exit
fi

############
#  OpenMP  #
############
 
TH=$2
 
# By default, without OpenMP
if [ -z $TH ]; then
 TH="none"
fi   
 
OPTM=""
 
if [ $TH = "none" ]; then
  if [ $FC = "ifort" ] || [ $FC = "mpiifort" ]; then
    #OPTM="-lmkl_sequential -lmkl_intel_lp64 -lmkl_core"
    OPTM="-lmkl_sequential"
  fi   
elif [ $TH = "omp" ]; then
  if [ $FC = "ifort" ] || [ $FC = "mpiifort" ]; then
    #OPTM="-qopenmp -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core"
    OPTM="-qopenmp -lmkl_intel_thread"
  else
    OPTM="-fopenmp"
  fi
else
  echo "Wrong OpenMP option ('omp' or 'none'). Exiting."
  exit
fi

OPT="$OPT $OPTM"

#################
#  Compilation  #
#################

echo "Starting the compilation process with $FC $OPT"
 
code=taurus_vap

# Creates a list of file to be compiled in the correct order
filelist="module_constants.xx module_mathmethods.xx MPI module_nucleus.xx \
          module_basis.xx module_hamiltonian.xx module_wavefunctions.xx   \
          module_fields.xx module_particlenumber.xx module_pairs.xx       \
          module_angularmomentum.xx module_multipoles.xx module_radius.xx \
          module_operators.xx module_projection.xx  module_constraints.xx \
          module_gradient.xx module_initialization.xx subroutines_pfaffian.yy"

# Remove the module parallelization when not doing MPI         
if [ $FC = "mpiifort" ] || [ $FC = "mpif90" ]; then
  filelist=$(echo "$filelist" | sed "s/MPI/module_parallelization.xx/g") 
else 
  filelist=$(echo "$filelist" | sed "s/MPI//g") 
fi
        
# The final list of .f90, .f and .o files
filef90=$(echo "$filelist" | sed "s/.xx/.f90/g" | sed "s/.yy/.f/g") 
fileo=$(echo "$filelist" | sed "s/.xx/.o/g" | sed "s/.yy/.o/g") 
filemod=$(echo "$filelist" | sed "s/.xx/.mod/g" | sed "s/module//g" \
                           | sed "s/\_//g")
         
# Creates wrk directory
wrkc=0
if [ ! -d $wrkdir ]; then 
  mkdir $wrkdir; echo "directory '$wrkdir' created"
  wrkc=1
fi

# Copy the files and removes mpi flag if necessary
for file in $filef90
do 
  cp $srcdir/$file $wrkdir/
  if [ $FC = "mpiifort" ] || [ $FC = "mpif90" ]; then 
    sed "s/\!cmpi /     /g" $wrkdir/$file > $wrkdir/tmp.f90
    mv $wrkdir/tmp.f90 $wrkdir/$file 
  fi
done

cp $srcdir/${code}.f90 $wrkdir/
if [ $FC = "mpiifort" ] || [ $FC = "mpif90" ]; then 
  sed "s/\!cmpi /     /g" $wrkdir/${code}.f90 > $wrkdir/tmp.f90 
  mv $wrkdir/tmp.f90 $wrkdir/${code}.f90
fi

echo "source files copied"

# Changes directory and performs the compilation
cd $wrkdir
         
for file in $filef90
do 
  echo "compiling ${file}"
  $FC $OPT -c $file
done

echo "compiling ${code}.f90"
$FC $OPT -o ${code}.exe ${code}.f90 $fileo $LIB

# Creates exe directory and move the exe file
if [ ! -d $exedir ]; then 
  mkdir $exedir; echo "directory '$exedir' created"
fi

if [ -f ${code}.exe ]; then mv ${code}.exe $exedir/; fi

##############
#  Clean up  #
##############

echo "cleaning up" 

cd $heredir

# Removes the wrkdir if not existing prior to the compilation
if [ $wrkc = 1 ]; then
 rm -rf $wrkdir
 echo "directory '$wrkdir' deleted"
else 
  for file in $filef90 $fileo $filemod
  do 
    rm -f $wrkdir/$file
  done
fi

# Final check to see if the exe was produced and move to exedir
if [ -f $exedir/${code}.exe ]; then 
  echo "compilation successful."
else
  echo "compilation failed."
fi
