# PhaseI Final Sar Pipeline

This is a effort to outline the final SAR 3D pipeline used for the TelEOSARus project along
with the scripts and Matlab files used.

## Outline of pipeline

This is the steps taken to generate the 3D SAR models from looking at the submission scripts

1. Spotlighting
 - Submission script: `run_fp###.sh` (pass)
 - Matlab Files: 
   - `Batch_Gotcha2006.m`
   - `SpotlightBasic.m`
 - Input Files: Original Gotcha data
 - Output Files: `blding620_FP###_0##.mat` (pass, number)

2. Backprojection
 - Submission script: `run_2d_fp###.sh`
 - Matlab Files: 
   - `Batch_ProceessData_Gotcha_phaseHistoryCorrectionJul12.m`
   - `formImageSimulationNew.m`
 - Input Files: `blding620_FP###_0##.mat`
 - Output Files: `full_2d_###_#.mat` (aperature, pass_idx)

3. Joint Sparse Recovery
 - Submission script: `run_full_3d_###.sh`
 - Matlab Files: 
   - `Batch_JointSparseRecovery_3D_script_Jul12.m`
   - `JointSparseRecovery_3D.m`
 - Input Files: `full_2d_###_#.mat` (aperature, pass_idx)
 - Output Files: `Results_3D_###.mat` (aperture)

4. Integrate over apertures 
 - Submission script: `run_integrate.sh`
 - Matlab Files: `image3d_integrate.m`
 - Input Files: `Results_3D_###.mat`
 - Output Files : 'resultsCombine.mat'

5. Create the point-cloud file in '.ply' format
 - Matlab Files: `pointCLoud_creator_parkingLot.m`
 - Input Files: `resultsCombine.mat`
 - Output Files : 'parkingLot_full_L1'
  
In case of the GOTCHA parking lot, skip steps 1 and 2 and 
proceed to Joint Sparse Recovery
 - Matlab Files: 
   - `JointSparseRecovery_3D_script.m`
   - `JointSparseRecovery_3D_L1.m`
 - Input Files: `saveData_#.mat` (pass_idx)
 - Output Files: `Results_3D_###.mat` (aperture)


