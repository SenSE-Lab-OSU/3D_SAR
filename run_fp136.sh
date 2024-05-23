#!/bin/bash

#PBS -l select=1:ncpus=48:mpiprocs=1
#PBS -l walltime=004:00:00
#PBS -q standard
#PBS -A WPrAFSNW27526A21
#PBS -m be
#PBS -M chet.nieter@kitware.com

module load matlab

for i in {0..150}; do
  matlab -nodisplay -nodesktop -r "pass=136;num=$i;run /p/home/cnieter/sar_image_code/Batch_Gotcha2006.m"
done

