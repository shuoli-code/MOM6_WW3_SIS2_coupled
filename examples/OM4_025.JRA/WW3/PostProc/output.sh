#!/bin/bash

#PBS -l ncpus=1
#PBS -l mem=192GB
#PBS -l jobfs=400GB
#PBS -q normal
#PBS -P groupname
#PBS -l walltime=1:00:00
#PBS -l storage=gdata/groupname+scratch/groupname
#PBS -l wd

module load openmpi
mpirun -n ${PBS_NCPUS} /build/intel/wave_ice_ocean/ww3_ounf/ww3_ounf
   