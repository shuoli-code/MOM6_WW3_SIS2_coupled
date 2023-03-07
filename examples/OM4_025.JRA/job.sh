#!/bin/bash
#PBS -l ncpus=240
#PBS -l mem=950GB
#PBS -l jobfs=2000GB
#PBS -q normal
#PBS -P groupname
#PBS -l walltime=8:00:00
#PBS -l storage=gdata/groupname+scratch/groupname
#PBS -l wd

module load openmpi
mpirun -n ${PBS_NCPUS} /MOM6_WW3_SIS2_coupled/build/intel/wave_ice_ocean/repro/MOM6

