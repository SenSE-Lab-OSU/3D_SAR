clc;
clear;
close all;
addpath('C:\Gotcha\spgl1-2.0');
addpath('C:\Gotcha\NUFFT_code');
cspeed = 299792458;

azSpan = 5; % degrees
numRangeBins = 1024;
numHeightBins = 256;
M_x=numRangeBins; M_y=numRangeBins; M_z=numHeightBins;
N  = [M_x M_y,M_z];
D_x = 100;
D_y = 100;
D_z = 120;

Res_x=D_x/M_x; Res_y=D_y/M_y; Res_z=D_z/M_z; 

lsize=D_x;

xImage=(-M_x/2+1:M_x/2)*Res_x;
yImage=(-M_y/2+1:M_y/2)*Res_y;

for azCenter =1:359
    
    for i=1:8
        passes = load(sprintf('saveData_%d',i));
        idxNeg =find(double(passes.th) < 0);
        passes.th(idxNeg)= 360 + passes.th(idxNeg);
        maxAz = mod(azCenter + azSpan/2,360);
        minAz = azCenter - azSpan/2;
        
        if (minAz < 0) || (minAz + azSpan >360)
            
            minAz = mod(360+minAz,360);
            
            idxAz1 =  find(double(passes.th) >= minAz);
            azimuthVals1{i} = double(passes.th(idxAz1));
            idxAz2 =  find(double(passes.th) < maxAz);
            azimuthVals1{i} = [azimuthVals1{i} double(passes.th(idxAz2))];
            idxAz = [idxAz1 idxAz2];
        else
            idxAz =  find((double(passes.th) < maxAz) & (double(passes.th) >= minAz));
            [ azimuthVals1{i},idxSorted] = sort(double(passes.th(idxAz)));
            idxAz = idxAz(idxSorted);
        end
        
        elev1{i} = double(passes.phi(idxAz));
        R0{i} = double(passes.r0(idxAz));
        f1{i}=repmat(double(passes.freq(:,1)),1,length(idxAz));
        numFreqs(i) = length(f1{i}(:,1));
        numPulses(i)= length(idxAz);
        p1 = (double(passes.fp(:,idxAz)));
        phCorrect{i} = double(passes.ph_correct(idxAz));
        rCorrect{i} = double(passes.r_correct(idxAz));
        curr_freq=2*pi*double(passes.freq(:,1));
        %ramp in seconds
        ramp=(2/cspeed)*rCorrect{i};
        %radial correction
        p1=p1.*exp(-1j*repmat(curr_freq,1,length(idxAz)).*repmat(ramp,length( f1{i}(:,1)),1));
        %phase correction
        p1=p1.*exp(1j*repmat(phCorrect{i},length( f1{i}(:,1)),1));
     
        
        AntX{i} = double(passes.x(idxAz));
        AntY{i} = double(passes.y(idxAz));
        AntZ{i} = double(passes.z(idxAz));
        
        taper_flag=0;
        data_pass(i) = formImageSimulationNew(elev1{i},numRangeBins,azimuthVals1{i},...
            p1,f1{i}(:,1),lsize,taper_flag,AntX{i},AntY{i},AntZ{i},...
            R0{i},xImage,yImage);
        im_final{i} = data_pass(i).im_final;
      %  data_pass(i).im_final =fliplr(flipud( data_pass(i).im_final));

        
    end
    clear passes;
    clear p1;
    clear idxNeg;
    clear phCorrect rCorrect;
    clear AntX AntY AntZ;
    
    xImage=(-M_x/2+1:M_x/2)*Res_x;
    yImage=(-M_y/2+1:M_y/2)*Res_y;

   JointSparseRecovery_3D_L1(im_final,xImage,yImage,f1,azimuthVals1,elev1);

end
