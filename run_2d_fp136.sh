#!/bin/bash

#PBS -l select=1:ncpus=48:mpiprocs=1
#PBS -l walltime=004:00:00
#PBS -q standard
#PBS -A WPrAFSNW27526A21
#PBS -m be
#PBS -M chet.nieter@kitware.com

module load matlab

matlab -nodisplay -nodesktop -r "pass=136;run /p/home/cnieter/sar_full_model_code/Batch_ProceessData_Gotcha_phaseHistoryCorrectionJul12.m"

