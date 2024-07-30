%% Prominent point processing method using the top-hat reflector
% Inputs 
% 1. topHat_coordinate_Center_x,topHat_coordinate_Center_y,topHat_coordinate_Center_z
% - The coordinates of the center of the top-hat reflector in the scene.
% 2. tophHat_radius - The radius of the top-hat reflector. The GOTCHA dataset
% has 1 m radius reflector. 
% 3. p1Input - Phase history measurements of dimension 
%               (Number of Frequencies X number of aperture points )
% 4. R0 - Row vector of the distances between scene center and the radar
% platform
% 5. f1 - frequencies used in each pulse  of size number of frequencies X
% number of pulses
% 6. azimuths - Azimuth angles used in the aperture.
% 7. AntX,AntY,AntZ - The coordinates of the aperture points. 

% Outputs
% 1. p1Final - Calibrated phase-history measurements
% 2. rCorrectEst - range correction  factor
% 3. phaseEst - phase correction factor

function [p1Final,rCorrectEst,phaseEst]=PPP_autofocus_twoPass(topHat_coordinate_Center_x,...
    topHat_coordinate_Center_y,topHat_coordinate_Center_z,tophHat_radius,...
    p1Input,R0,f1,azimuths,AntX,AntY,AntZ)

p1 =p1Input;
cspeed = 299792458;
bandwidth  = f1(end,1)-f1(1,1);
oversamplingFactor = 16;
range_res = cspeed/(2*bandwidth*oversamplingFactor);
numFreqs = length(f1(:,1))*oversamplingFactor;
searchWindow = 2;
delR =(-numFreqs/2:numFreqs/2-1)*range_res;

numPulses= length(azimuths);

% The coordinates of top-hat dominant scattering center move along with the 
% radar platform.
topHat_coordinate_x =  + tophHat_radius*cosd(azimuths);
topHat_coordinate_y = + tophHat_radius*sind(azimuths);
topHat_coordinate_z=zeros(size(topHat_coordinate_y));

% Differential range of the top-hat coordinates with respect to scene
% center
diffRangeTrue = (sqrt(sum(([ topHat_coordinate_x.' topHat_coordinate_y.' topHat_coordinate_z.'] - ...
    [AntX.' AntY.' AntZ.']).^2,2))-R0.');
% Differential range of the top-hat center with respect to scene
% center
diffRangeCenter = (sqrt(sum(([ topHat_coordinate_Center_x.' topHat_coordinate_Center_y.' topHat_coordinate_Center_z.'] - ...
    [AntX.' AntY.' AntZ.']).^2,2))-R0.');
RangeCenter = (sqrt(sum(([ 0 0 0] - ...
    [AntX.' AntY.' AntZ.']).^2,2)));
RangetopHat = (sqrt(sum(([ topHat_coordinate_x.' topHat_coordinate_y.' topHat_coordinate_z.'] - ...
    [AntX.' AntY.' AntZ.']).^2,2)));

%Shifting the scene center to the top-hat center
curr_freq=2*pi*f1;
rampCenterShift=(2/cspeed)*diffRangeCenter.';
p1 = p1.*exp(1j*curr_freq.*repmat(rampCenterShift,length( f1(:,1)),1));
% p1Calib = p1Calib.*exp(1j*repmat(curr_freq,1,numPulses).*repmat(rampCenterShift,length( f1),1));

%Computing the range-profile
pp = ifftshift(ifft(p1,numFreqs),1);

% Correcting the range error based on the expected location of the peak
% scattering center corresponding to the top-hat.


for j=1:size(pp,2)
    idxR =find((delR < diffRangeTrue(j) +searchWindow) & (delR > diffRangeTrue(j) -searchWindow) );
    ppSub = pp(idxR,j);  
    ppSub= [min(ppSub); ppSub; min(ppSub)];

    delRSub = delR(idxR).';
    eee=delRSub(1);
    fff= delRSub(end);
    delRSub = [eee-0.01; delRSub;fff+0.01 ];
    [ii,kk]=findpeaks(20*log10(abs(ppSub)),'MinPeakHeight',max(20*log10(abs(ppSub)))-2 );
    
    pkphase(j)=angle(max(ppSub));
    [~,jj] = min(abs(delRSub(kk)-diffRangeTrue(j)));
    
    diffRangeEst(j)= delRSub(kk(jj));
    phaseAdd(j) = mod(4*pi*mean(f1(:,j))/cspeed*(RangeCenter(j)-RangetopHat(j)),2*pi)-pi;
end

% Range correction 
rCorrectEst = -(diffRangeEst.'-diffRangeTrue);
rCorrectMedian = median(rCorrectEst);
ramp=(2/cspeed)*rCorrectMedian*ones(size(rCorrectEst.'));
p1=p1.*exp(-1j*curr_freq.*repmat(ramp,length( f1(:,1)),1));
pp = ifftshift(ifft(p1,numFreqs),1);

searchWindowNew = 0.25;

for j=1:size(pp,2)
    idxR =find((delR < diffRangeTrue(j) +searchWindowNew) & (delR > diffRangeTrue(j) -searchWindowNew) );
    ppSub = pp(idxR,j);  
    ppSub= [min(ppSub); ppSub; min(ppSub)];

    delRSub = delR(idxR).';
    eee=delRSub(1);
    fff= delRSub(end);
    delRSub = [eee-0.01; delRSub;fff+0.01 ];
    [ii,kk]=findpeaks(20*log10(abs(ppSub)),'MinPeakHeight',max(20*log10(abs(ppSub)))-2 );
    
    pkphase(j)=angle(max(ppSub));
    [~,jj] = min(abs(delRSub(kk)-diffRangeTrue(j)));
    
    diffRangeEst(j)= delRSub(kk(jj));
end

rCorrectEst = (diffRangeEst.'-diffRangeTrue);

ramp=(2/cspeed)*rCorrectEst.';
% Apply the range-correction
p1=p1.*exp(-1j*curr_freq.*repmat(ramp,length( f1(:,1)),1));



pp = ifftshift(ifft(p1,numFreqs),1);

% Estimate the required phase correction 
for j=1:size(pp,2)
    idxR =find((delR < diffRangeTrue(j) +searchWindow) & (delR > diffRangeTrue(j) -searchWindow) );
    ppSub = pp(idxR,j);
    ppSub= [min(ppSub); ppSub; min(ppSub)];
    delRSub = delR(idxR).';
    eee=delRSub(1);
    fff= delRSub(end);
    delRSub = [eee-0.01; delRSub;fff+0.01 ];
    [ii,kk]=findpeaks(20*log10(abs(ppSub)),'MinPeakHeight',max(20*log10(abs(ppSub)))-2 );
    
    phaseEst(j) = -angle(max(ppSub)) + phaseAdd(j);
end

rCorrectEst = rCorrectMedian + rCorrectEst;
%ramp in seconds
ramp=(2/cspeed)*rCorrectEst.';
%radial correction
p1Final=(p1Input).*exp(-1j*curr_freq.*repmat(ramp,length( f1(:,1)),1));
%phase correction
p1Final=p1Final.*exp(1j*repmat(phaseEst,length( f1(:,1)),1));
end