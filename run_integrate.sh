#!/bin/bash

#PBS -l select=1:ncpus=48:mpiprocs=1
#PBS -l walltime=001:00:00
#PBS -q debug
#PBS -A WPrAFSNW27526A21
#PBS -m be
#PBS -M chet.nieter@kitware.com

module load matlab

matlab -nodisplay -nodesktop -r "run /p/home/cnieter/sar_full_model_code/image3d_integrate.m"

