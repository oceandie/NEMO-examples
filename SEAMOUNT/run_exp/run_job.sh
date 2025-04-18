#!/bin/bash -l
#
#PBS -P climate
#PBS -N SEAMOUNT
#PBS -q normal
#PBS -l walltime=1800
#PBS -l select=1

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

ulimit -c unlimited
ulimit -s unlimited

NEMO_N=4
XIOS_N=1

echo " mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./nemo : -n $XIOS_N -d 1 ./xios"
mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./nemo : -n $XIOS_N -d 1 ./xios

