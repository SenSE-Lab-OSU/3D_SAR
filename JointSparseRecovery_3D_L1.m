%% Generates the 3D representation and saves the result in a file
% inputs
% im_final -  number of passes x 1 cell array containing the images from
%               all the passes
% xImage - image coordinates in x direction
% yImage - image coordinates in y direction
% f1 -  number of passes x 1 cell array containing Frequencies 
%               used in each pass
% azimuthVals -  number of passes x 1 cell array containing  
%               azimuth angles of the flight path in each pass
% elev1 -  number of passes x 1 cell array containing  
%               elevation angles of the flight path in each pass


function JointSparseRecovery_3D_L1(im_final,xImage,yImage,f1,azimuthVals,elev1,idxImage)

cspeed = 299792458;

snrReconstruction= 5; %(dB)
% SNR threhsold for cleaning up the clutter
snrThreshold = 80; %(dB)

numPasses = length(im_final);
azCenter = round(median(azimuthVals{1}));

numRangeBins = length(xImage);
% This refers to the number of voxels in z-dimension and X-Y plane 
numHeightBins = 256; %z
numRangeBinsVoxel = 1024; %X-Y
M_x=numRangeBinsVoxel;
M_y=numRangeBinsVoxel;
M_z=numHeightBins;

D_x = abs(xImage(end) - xImage(1));
D_y = abs(yImage(end) - yImage(1));
D_z = 80; % height of scene [-40m,40m], change this according to scene size
shiftZ=D_z/2-10; % shift in z dimension in meters
Res_x=D_x/numRangeBins; Res_y=D_y/numRangeBins;
Res_xVoxel=D_x/M_x; Res_yVoxel=D_y/M_y; Res_z=D_z/M_z; 


xImageVoxel=(-M_x/2+1:M_x/2)*Res_xVoxel;
yImageVoxel=(-M_y/2+1:M_y/2)*Res_yVoxel;



zImage=(-M_z/2+1:M_z/2)*Res_z;

%[zGrid,~,~] = ndgrid(zImage,yImageVoxel,xImageVoxel);
% groups1=zeros(size(zGrid));
% 
% for i=1:size(zGrid,1)
%     groups1(i,:,:) = reshape(single(1:M_x*M_y).',1,M_y,M_x);
% end
% groups = single(groups1(:));
% clear zGrid;
% clear groups1;


for i=1:numPasses
    numFreqs(i) = length(f1{1}(:,1));
    numPulses(i)= length(azimuthVals{i});    
    im_final{i} =fliplr(flipud(im_final{i}));
    meanF(i)= mean(f1{i}(:,1));
    meanElev(i) = mean(elev1{i});
end

[xGrid1,yGrid1] = meshgrid(xImage,yImage);

fixedCenterFreq =  mean(meanF);
fixedElev =  mean(meanElev);

phaseRampCenter_y = 2*cosd(fixedElev)*sind(azCenter)*fixedCenterFreq/cspeed;
phaseRampCenter_x = 2*cosd(fixedElev)*cosd(azCenter)*fixedCenterFreq/cspeed;
phaseRampCenter_z = 2*sind(fixedElev)*fixedCenterFreq/cspeed;

kx_grid=linspace(-1/2/Res_x, 1/2/Res_x, numRangeBins+1); kx_grid(end)=[];
ky_grid=linspace(-1/2/Res_y, 1/2/Res_y, numRangeBins+1); ky_grid(end)=[];

kx_grid_voxel=linspace(-1/2/Res_xVoxel, 1/2/Res_xVoxel, M_x+1); kx_grid(end)=[];
ky_grid_voxel=linspace(-1/2/Res_yVoxel, 1/2/Res_yVoxel, M_y+1); ky_grid(end)=[];

kz_grid=linspace(-1/2/Res_z, 1/2/Res_z, M_z+1); kz_grid(end)=[];

k_x_total=[];
k_y_total=[];
k_z_total=[];
phTotal=[];
samplesIndexPass=1;
for i=1:numPasses
    im_final{i} = im_final{i}...
        .*exp(-1i*2*pi*(phaseRampCenter_y*yGrid1+phaseRampCenter_x*xGrid1));
    im_final{i} = im_final{i}.';
    samplesIndexPass(i+1)=length(f1{i}(:,1))*length(azimuthVals{i});
    frep = f1{i};%repmat(f1{i},1,length(azimuthVals{i}));
    azrep = repmat(((azimuthVals{i})),length(f1{i}(:,1)),1);
    elrep = repmat(((elev1{i})),length(f1{i}(:,1)),1);
    
    k_y=2/cspeed*sind(azrep).*(frep).*cosd(elrep);
    k_y=k_y(:);
    k_y_c = phaseRampCenter_y;
    k_y = k_y -k_y_c  ;
    
    k_x=2/cspeed*cosd(azrep).*(frep).*cosd(elrep);
    k_x=k_x(:);
    k_x_c =phaseRampCenter_x;
    k_x = k_x - k_x_c;

    k_z=2/cspeed*(frep).*sind(elrep);
    k_z_c =phaseRampCenter_z;
    k_z=k_z(:) -k_z_c;
    im_final{i}  = fliplr(flipud(im_final{i}));
    [y,k1] = iFGG_2d_type2mod(im_final{i},[k_x k_y],12,kx_grid,ky_grid);
    y=reshape(y,size(f1{i}));
    y = exp(-1i*2*pi*reshape(k_z,size(f1{i}))*shiftZ).*y;
    y=y(:);
    y=y/norm(y);
    phTotal=[phTotal;y];
    k_x_total=[k_x_total;k_x];
    k_y_total=[k_y_total;k_y];
    k_z_total=[k_z_total;k_z];
    k_x_1{i} = k_x;
    k_y_1{i} = k_y;
    k_z_1{i} = k_z;
    
end
clear azrep;
clear elrep;
clear frep;
clear k_x;
clear k_y;
clear k_z;
clear data_pass;
%% Solving the regularized image-reconstruction for all the passes jointly as a 3D problem
samplesIndexPass=cumsum(samplesIndexPass)-1;
samplesIndexPass(1)=0;

phTotal=numPasses*phTotal/norm(phTotal);

sigma_n = numPasses*10^(-snrReconstruction/10);

A1=@(x,mode)sar_operator_nufft_3d(x,mode,k_x_total,k_y_total,k_z_total,...
    kx_grid_voxel,ky_grid_voxel,kz_grid,M_x,M_y,M_z);

options = spgSetParms('isComplex',1,'verbosity',1,'optTol',7e-3);
X2=  spg_bpdn(A1,phTotal, sigma_n, options ); % Changed this to solve sparsity promoting regularizer
clear groups;

[zGrid,yGrid,xGrid] = ndgrid(zImage,yImageVoxel,xImageVoxel);
p_x=[];
p_y=[];
p_z=[];
c_1=[];

%% Performing the 3D reconstruction by reshaping the obtained result
xx = reshape(X2,M_z,M_y,M_x);
xx= permute(xx,[1,3,2]);
 


[aa,bb,cc]=ind2sub(size(xx),find(20*log10(abs(xx)) > max(20*log10(abs(xx(:)))) - snrThreshold));
sc_points_layOver(1,:) = arrayfun(@(i)xGrid(aa(i),bb(i),cc(i)),(1:length(aa)).');
sc_points_layOver(2,:) = arrayfun(@(i)yGrid(aa(i),bb(i),cc(i)),(1:length(aa)).');
sc_points_layOver(3,:)=  shiftZ + arrayfun(@(i)zGrid(aa(i),bb(i),cc(i)),(1:length(aa)).');

amps= 20*log10(abs(arrayfun(@(i)xx(aa(i),bb(i),cc(i)),(1:length(aa)).')));
viewAngle = azCenter;

% figure; scatter3(sc_points_layOver(1,:) ,sc_points_layOver(2,:),...
%     sc_points_layOver(3,:),40,((amps)),'filled');
% colormap(flipud(SAR_cmap)); colorbar;

save(sprintf('/p/work1/cnieter/sar_images/3d_results/Results_3D_%d',idxImage),'fixedElev','fixedCenterFreq',...
    'meanElev','meanF','azCenter','amps','sc_points_layOver','viewAngle','-v7.3');

end
