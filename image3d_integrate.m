clc;
clear;
close all;


d = dir('/p/work1/cnieter/sar_images/3d_results/Results_3D_*.mat');
snrThreshold = 50; % dB below max value that is diplayed

numImages = length(d);
load(sprintf('/p/work1/cnieter/sar_images/3d_results/%s', d(1).name));

pointsTotal = sc_points_layOver.';
ampsTotal = amps;

for idxImages = 2:numImages
    load(sprintf('/p/work1/cnieter/sar_images/3d_results/%s', d(idxImages).name));
    % Find the intersection of the points from all previous viewws with the
    % new view
    [pointsTemp,idxTotal,idxNew] = intersect(pointsTotal,sc_points_layOver.','rows'); 
    
    % combine all the amplitudes of the common points
    ampsCombine = 20*log10(10.^(ampsTotal(idxTotal)/20)+10.^(amps(idxNew)/20));
    % replace the resulting superposition at the original locations
    ampsTotal(idxTotal) = ampsCombine;
    
    % Find the unique points in the new view from all the previous views
    [pointsTempDiff,idxNewDiff] = setdiff(sc_points_layOver.',pointsTotal,'rows');
    
    % append the newly found points
    pointsTotal = [pointsTotal;pointsTempDiff];
    ampsTotal = [ampsTotal;amps(idxNewDiff)];
    fprintf(sprintf('Finished file %s', d(idxImages).name))
end

% find the points that have amplitudes above a threshold
ampsTotalFinal = ampsTotal(ampsTotal >max(ampsTotal) -snrThreshold);
pointsTotalFinal = pointsTotal(ampsTotal >max(ampsTotal) -snrThreshold,:);

save('/p/work1/cnieter/sar_images/3d_results/resultsCombine','ampsTotalFinal','pointsTotalFinal','-v7.3');
