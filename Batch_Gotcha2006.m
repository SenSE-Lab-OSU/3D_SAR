
% Input parameters
% num = 0;
% pass = 136;

aperture = 8192;
pulsesPerFile = 27092;

% Region of interest
locX= 200;
locY= -290;
locZ= 11.41;
vX=250;
vY=250;

% Define the path to the base directory of the dataset
% these parameters/variables are all dataset dependent)
dataPath = sprintf('/p/work1/cnieter/gotcha/2006/MOCOMP_GD/FP0%d/MCPH/HH_C01', pass);
if pass < 140
    filebase = sprintf('C1SAP1C0100%d', pass + 381); % Base filename 
else
    filebase = sprintf('C1SAP1C0100%d', pass + 378); % Base filename 
end
startPulse = num*aperture; % Starting pulse

% Path to save images
resultPath = sprintf('/p/work1/cnieter/sar_images/FP0%d/', pass);

% initiate variables
data.minF = [];
data.AntX = [];
data.AntY = [];
data.AntZ = [];
data.R0 = [];
data.phdata = [];
data.Np = aperture;

% File parameters
fileParams = {};

% Find first file to process 
ff.start = mod(startPulse, pulsesPerFile);
ff.fileIdx = idivide(int64(startPulse), int64(pulsesPerFile));
fileParams{1} = ff; 

% Check if second file is needed
if fileParams{1}.start + aperture > pulsesPerFile
    sf.start = 0;
    sf.Np = mod(startPulse + aperture, pulsesPerFile);
    sf.fileIdx = idivide(int64(startPulse + aperture), int64(pulsesPerFile));
    fileParams{2} = sf; 
    fileParams{1}.Np = pulsesPerFile - fileParams{1}.start;
else
    fileParams{1}.Np = aperture;
end

for i=1:length(fileParams)

    fnamePHdata = sprintf('%s%s%s_%03d',dataPath,filesep,filebase,fileParams{i}.fileIdx);
    fnamePHXheader = sprintf('%s.phxhdr',fnamePHdata);
    fnamePAUX = sprintf('%s%s%s.paux',dataPath,filesep,filebase);

    % Read in the phoenix header to gather needed parameters
    fid = fopen(fnamePHXheader,'r');
    while 1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        I = find(tline=='=');
        switch tline(1:(I-1))
            % Number of samples per pulse
            case('NumberOfSamples')
                data.K = str2double(tline((I+1):end));

                % The bandwidth of each pulse
            case('SensorBandwidth')
                dtemp = tline((I+1):end);
                J = find(dtemp==' ');
                BW = str2double(dtemp(1:(J(2)-1))) * 1e6;
        end
    end
    fclose(fid);

    % Calculate frequency step size
    data.deltaF = BW/data.K;

    % Save data from header file
    data.J = J;
    data.BW = BW;

    if length(data.minF) == 0
        % Read in PAUX data
        paux = readPostDCPHPaux(fnamePAUX, data.Np, startPulse);

        % Minimum frequency for each pulse
        data.minF = [data.minF; paux(:,11)];

        % (x,y,z) coordinates of the antenna at each pulse (m)
        data.AntX = [data.AntX; paux(:,17)];
        data.AntY = [data.AntY; paux(:,18)];
        data.AntZ = [data.AntZ; paux(:,19)];

        % Distance to scene center for each pulse (m)
        data.R0 = [data.R0; paux(:,20)];
    end

    % Read in Phase History Data
    data.phdata = [data.phdata, readPostDCPH(fnamePHdata, data.K, fileParams{i}.Np, fileParams{i}.start).'];

    fprintf('\nfinished file: %s\n', fnamePHdata);
end

% Clip off zeros in phase data
minind=zeros(data.Np,1);maxind=zeros(data.Np,1);

for ii=1:data.Np
    minind(ii)=min(find(abs(data.phdata(:,ii))>eps));
    maxind(ii)=max(find(abs(data.phdata(:,ii))>eps));
end;

data.K=round (mean(maxind-minind)+1);

phdata_old=data.phdata;
data.phdata=zeros(data.K, data.Np);

for ii=1:data.Np
    data.phdata(:,ii)=phdata_old(minind(ii):minind(ii)+data.K-1,ii);
end

orig_size = [data.K,data.Np];
data.deltaF=data.BW/data.K;
% Spotlight data
fprintf(1,'Spotlighting data...\n');
Center = [locX locY locZ];
View = [vX vY];
data = SpotlightBasic(data,Center,View);

data.orig_size = orig_size;

% Save the data
resultName = sprintf('blding620_FP%03d_%03d', pass, num);
save(strcat(resultPath,resultName),'data');
fprintf('Saved file %s/n', strcat(resultPath,resultName))

fprintf('Finished with section %d\n',num);

exit
