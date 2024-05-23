%clc;
%clear;
%close all;

pathnameSPGl1 = '/p/home/cnieter/sar_full_model_code/spgl1-2.0/';
pathnameNUFFT = '/p/home/cnieter/sar_full_model_code/NUFFT_code/';

addpath(pathnameSPGl1);
addpath(pathnameNUFFT);

% This is the number of elevation passes. Gotcha dataset has 8 passes
numPasses = 8; 
resultsName = '/p/work1/cnieter/sar_images/2d_image_files/full_2d';

for idxImages= start:stop
    for idxPass =1:numPasses
        d1(idxPass) = load(sprintf('%s_%d_%d',resultsName,idxImages,idxPass));
        xImage = d1(1).xImage;
        yImage = d1(1).yImage;
        im_final{idxPass} = d1(idxPass).im_final;
        f1{idxPass} = d1(idxPass).freqFinal;
        azimuthVals{idxPass} = d1(idxPass).thFinal;
        elev1{idxPass} = d1(idxPass).phiFinal;
    end
    JointSparseRecovery_3D_L1(im_final,xImage,yImage,f1,azimuthVals,elev1,idxImages);
end

exit

