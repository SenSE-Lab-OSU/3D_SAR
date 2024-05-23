
% clear;
% close all;

% Assuming each pass in a different folder
folders = sprintf('/p/work1/cnieter/sar_images/FP0%d',pass)
% Get all the mat-files names
d1=dir(sprintf('%s%s*.mat',folders,filesep));
% Extracting the first part from the file name ignoring the look-angle part
% extracting 'blding620_FP146' from 'blding620_FP146_001.mat'

tempStr= cellfun(@(x)x(1:end-8),{d1.name},'UniformOutput',false);
uniqueStr = unique(tempStr);

%% Top-hat details
topHat_coordinate_Center_x = 52.32;
topHat_coordinate_Center_y = 121.5 %122.39;
topHat_coordinate_Center_z = 0;

tophHat_radius=1;
AutoFocusEnable = 1;


%% scene-details
sceneExtent_x = 350;
sceneExtent_y = 350;

%% Image-details
NFFT = 131072;
x0 =0;  % center of image
y0= 0;   % center of image


Nx = 1024; % number of pixels in x direction
Ny = 1024; % number of pixels in y direction


Res_x = sceneExtent_x/Nx;
Res_y = sceneExtent_y/Ny;

x_min = -sceneExtent_x/2;
x_max = sceneExtent_x/2 ;

y_min = -sceneExtent_y/2;
y_max = sceneExtent_y/2 ;


xImage=x_min+Res_x:Res_x:x_max;
yImage=y_min+Res_y:Res_y:y_max;


%The images will be stored in a file named resultsName_passNuber.mat
resultsName = strcat(folders, '/full_2d');
 

%% Performing auto-focus on the phase-history measurements
for idxPass = 1:length(uniqueStr)

    phdata=[]; % phase-history data after auto-focus
    freq=[]; % frequencies used
    x=[]; % x-location of radar
    y=[]; % y-location of radar
    z=[]; % z-location of radar
    r0=[]; % scene center
    th=[]; % azimuth
    phi=[]; % elevation
    rCorrectEstSave = []; %range-correction factor
    phaseEstSave = []; % phase-correction factor
    phdataOriginal = []; % original phase-history meassurements
    
    % extract the indices of files from the current pass
    indexFiles =  find(contains({d1.name},uniqueStr{idxPass})==1);
    
    for idxFiles= 1:length(indexFiles)
        fprintf('prcessing image %d in pass %d\n',idxFiles,idxPass);
        
        load(sprintf('%s%s%s',d1(indexFiles(idxFiles)).folder,filesep,d1(indexFiles(idxFiles)).name));
        
        numFreqSamples = size(data.phdata,1);
        numPulses = size(data.phdata,2);
        
        if numFreqSamples < 1981
            data.deltaF=3.231472355518341e+05;
            data.phdata=[data.phdata ;zeros(1981-numFreqSamples,numPulses)];
        elseif numFreqSamples > 1981
                data.phdata=data.phdata(1:1981,:);
        end
        
        numFreqSamples =1981;
        data.K=1981;
        data.deltaF=  3.231472355518341e+05;
        
        freqTemp =double(data.minF) + data.deltaF*(0:numFreqSamples-1).';
        freq = [freq freqTemp];
        thTemp=double(data.azim);
       % thTemp(thTemp <0) =thTemp(thTemp <0) + 2*pi;  % making angles between 0 to 2*pi
        
        x = [x double(data.AntX)];
        y = [y double(data.AntY)];
        z = [z double(data.AntZ)];
        r0 = [r0 double(data.R0)];
        th = [th thTemp]; % assuming azimuth angles are in radians
        phi = [phi rad2deg(double(data.elev))]; % assuming elevation angle are in radians
        
        %% performing auto-focus using prominent point method if flag is enabled
        if (AutoFocusEnable == 1)
            [p1Final,rCorrectEst,phaseEst]=PPP_autofocus_twoPass(topHat_coordinate_Center_x,...
                topHat_coordinate_Center_y,topHat_coordinate_Center_z,tophHat_radius,...
                double(data.phdata), double(data.R0),freqTemp,rad2deg(thTemp),double(data.AntX),...
                double(data.AntY),double(data.AntZ));
            plot(rCorrectEst);
            phaseEstSave = [phaseEstSave phaseEst];
            rCorrectEstSave = [rCorrectEstSave rCorrectEst.'];
            phdata = [phdata p1Final];
        else
            phdata = [phdata double(data.phdata)];
        end
      
        phdataOriginal =  [phdataOriginal double(data.phdata)];
        if idxFiles == 1 % store the starting angle from file number 1
            startingTh =thTemp(1);
        end
        
    end



th = rad2deg(unwrap(th)); % unwrapping the angles to see if there is any mismatch

%     selecting starting angle  to starting angle + 360 degree in azimuth
 %check if it is decreasing  in azimuth, if not flip everything
 
 if th(1)< th(end)
     idx=length(th):-1:1;
     th = th(idx);
    freq = freq(:,idx);
    x = x(idx);
    y = y(idx);
    z = z(idx);
    r0 = r0(idx);
    phi = phi(idx);
    phdata = phdata(:,idx);
    phdataOriginal = phdataOriginal(:,idx);
 end
 
     
startingTh=th(1);
idx = find(th <= startingTh & th > startingTh - 360 );
th = th(idx);
freq = freq(:,idx);
x = x(idx);
y = y(idx);
z = z(idx);
r0 = r0(idx);
phi = phi(idx);
phdata = phdata(:,idx);
phdataOriginal = phdataOriginal(:,idx);

if (AutoFocusEnable == 1)
    phaseEstSave = phaseEstSave(idx);
    rCorrectEstSave = rCorrectEstSave(idx);
end
save(sprintf('%s/phaseHistory_correction_%d_pass',folders,idxPass),'phaseEstSave','rCorrectEstSave','-v7.3');

th = mod(th,360); % wrapping it back to 360

%% Constructing images using the auto-focussed data.
%creating 5-degreee images
numImages = 144; %  number of images needed
azSkip = 360/numImages;

azCenter = 0:azSkip:360-azSkip;

azSpan = 5; % azSpan degree images will be created
for idxImages =1:numImages
    azMin = mod(azCenter(idxImages) - azSpan/2,360);
    azMax = mod(azCenter(idxImages) + azSpan/2,360);
    idx1 = find(th < azMax);
    idx2 = find(th > azMin);
    if (azCenter(idxImages) - azSpan/2 < 0) ||  (azCenter(idxImages) + azSpan/2 >= 360)
        % handle the boundaries
        idxFinal= [idx2 idx1];
    else
        % the main part 
        idxFinal = intersect(idx1,idx2);
    end
    thImage = th(idxFinal);
    freqImage = freq(:,idxFinal);
    AntXImage = x(idxFinal);
    AntYImage = y(idxFinal);
    AntZImage = z(idxFinal);
    r0Image = r0(idxFinal);
    phiImage = phi(idxFinal);
    phdataImage = phdata(:,idxFinal);
    taper_flag =0; % windowing disabled
    
    % forming the image using Wx, Wy, Nfft, x_vec, y_vec from the input
    % data file and assuming it to be constant for all files in the pass.
    
    data_pass = formImageSimulationNew(phiImage,Nx,Ny, thImage,...
        phdataImage,freqImage(:,1),sceneExtent_x,sceneExtent_y,NFFT,taper_flag,...
        AntXImage,AntYImage,AntZImage,r0Image,xImage,yImage,x0,y0);
    im_final = data_pass.im_final;
    thFinal = thImage;
    freqFinal = freqImage;
    AntXFinal= AntXImage;
    AntYFinal= AntYImage;
    AntZFinal= AntZImage;
    r0Final = r0Image;
    phiFinal = phiImage;
    save(sprintf('%s_%d_%d',resultsName,idxImages,idxPass),'freqFinal','AntXFinal',...
        'AntYFinal','AntZFinal','r0Final','thFinal','phiFinal','im_final',...
        'xImage','yImage','-v7.3');
end

end

exit
