function data=SpotlightBasic(data,Center,View)
%******************************************************************************
% This function digitally spotlights phase history from a larger data
% source. The data is compressed using frequency (fast time) and pulse (slow
% time)decimation while preserving alias-free range and crossrange in the
% specified View. All sensor location and range information is referenced to
% the new scene center in the output. This function assumes that data.minF and
% data.deltaF are the same for all pulses.
%
% Inputs...
% data.phdata: phase history data, fast time in rows, slow time in columns
% data.deltaF: step size of frequency data (Hz), assume all are the same
% data.minF: vector containing the start frequency (Hz) of each pulse
% data.AntX: x-position (m) of the sensor at each pulse
% data.AntY: y-position (m) of the sensor at each pulse
% data.AntZ: z-position (m) of the sensor at each pulse
% data.R0: range (m) from scene center to sensor for each pulse
% data.K:  the number of frequency bins per pulse
% data.Np: the number of pulses
% Center: spotlight center is [x,y,z] (m) with respect to original scene center
% View: spotlight shall enable alias-free imaging in a scene of size [x,y] (m)
%
% Outputs...
% data: same structural variable as input
% data.View: store the View input for informative purposes
%******************************************************************************

%set parameters
SPT=~isempty(ver('signal')); %detect signal processing toolbox
c=299792458; %speed of light in m/s
derate=1.05; %make spot slightly larger than necessary
data.View=View; %store the view information

%some functions require doubles but data may be in single
data.phdata=double(data.phdata); data.deltaF=double(data.deltaF);
data.minF=double(data.minF); data.AntX=double(data.AntX);
data.AntY=double(data.AntY); data.AntZ=double(data.AntZ);
data.R0=double(data.R0); data.K=double(data.K); data.Np=double(data.Np);

%******************************************************************************
% re-center and decimate in fast time
%******************************************************************************

%calculate the frequency decimation factor for alias-free range
deltaFspot=c/(2*derate*norm(data.View));
N=floor(deltaFspot/data.deltaF);

for n=1:data.Np %for each pulse
    
    %re-center
    freq=(data.minF(n)+data.deltaF*(0:(data.K-1))).';
    R0new=norm([data.AntX(n)-Center(1),data.AntY(n)-Center(2),data.AntZ(n)-Center(3)]);
    dR=data.R0(n)-R0new;
    data.phdata(:,n)=data.phdata(:,n).*exp(-1i*4*pi/c*dR*freq);
    data.R0(n)=R0new; %update range to new center
    data.AntX(n)=data.AntX(n)-Center(1);
    data.AntY(n)=data.AntY(n)-Center(2);
    data.AntZ(n)=data.AntZ(n)-Center(3);
    
    %decimate frequency bins
    freq=decimate_(freq,N,SPT);
    data.minF(n)=freq(1);
    
    %decimate PH in frequency
    phdata_(:,n)= decimate_(real(data.phdata(:,n)),N,SPT)+...
        1i*decimate_(imag(data.phdata(:,n)),N,SPT);
    
end

%store the decimated data
data.phdata=phdata_;
data.K=size(data.phdata,1);
data.deltaF=freq(2)-freq(1); %assume that deltaF is the same for all pulses

%******************************************************************************
% interpolate the PH and Ant information to uniform azimuth spacing
%******************************************************************************

%analyze the azimuth spacing of pulses
[Ant.az,Ant.el,Ant.R] = cart2sph(data.AntX,data.AntY,data.AntZ);
Ant.az=unwrap(Ant.az);
RPP=Ant.az(2:data.Np)-Ant.az(1:(data.Np-1)); %radians per pulse, directional
abs_RPP=abs(RPP); %magnitude of radians per pulse
[sort_RPP,I]=sort(abs_RPP);
im=I(5); RPPdata=sort_RPP(5); %ignor up to 4 anomalies

%interpolate the PH and Ant to uniform grid with min_RPP spacing
az_i=Ant.az(1):RPP(im):Ant.az(end);
el_i=interp1(Ant.az,Ant.el,az_i);
R_i=interp1(Ant.az,Ant.R,az_i);
data.Np=size(az_i,2);
PH_i=zeros(data.K,data.Np);
for k=1:data.K, PH_i(k,:)=interp1(Ant.az,data.phdata(k,:),az_i); end
data.phdata=PH_i;

%update the antenna and range using the interpolated locations
Ant.az=az_i; Ant.el=el_i; Ant.R=R_i;
[data.AntX,data.AntY,data.AntZ]=sph2cart(Ant.az,Ant.el,Ant.R);
data.R0=Ant.R;

%******************************************************************************
% decimate in slow time
%******************************************************************************

%calculate the PPR decimation factor for alias-free crossrange
fmax=data.minF(1)+(data.K-1)*data.deltaF; %find the maximum frequency bin
PPRspot=derate*2*norm(data.View)*fmax*cos(min(Ant.el))/c; %new PPR
PPRdata=1/RPPdata; %old PPR
M=floor(PPRdata/PPRspot);

%decimate the phase history in slow time
FilterScale=decimate_(ones(size(Ant.az)),M,SPT); %detect the filter scaling
for k=1:data.K %across each frequency bin
    phdataNew(k,:)=(decimate_(real(data.phdata(k,:)),M,SPT)+...
        1i*decimate_(imag(data.phdata(k,:)),M,SPT))./FilterScale;
end
%decimate the antenna location in slow time
AntNew.az=decimate_(Ant.az,M,SPT)./FilterScale;
AntNew.el=decimate_(Ant.el,M,SPT)./FilterScale;
AntNew.R=decimate_(Ant.R,M,SPT)./FilterScale;

data.azim = AntNew.az; 
data.elev = AntNew.el; 
data.R = AntNew.R; 

%update the output data structural variable
data.phdata=phdataNew;
data.Np=size(data.phdata,2);
[data.AntX,data.AntY,data.AntZ]=sph2cart(AntNew.az,AntNew.el,AntNew.R);
data.R0=AntNew.R;
data.minF=ones(1,data.Np)*data.minF(1,1);

end

%******************************************************************************
% Select between the decimate function in the Matlab Signal Processing
% Toolbox and the basic decimate_alt function that requires no toolbox.
%******************************************************************************
function Y=decimate_(X,M,SPT)
if M>1
    if SPT
        Y=decimate(X,M,'FIR'); %in Matlab's signal processing toolbox
    else
        Y=decimate_alt(X,M); %an alternate decimator
    end
else
    Y=X;
end

end

%******************************************************************************
% The decimate_alt function implements a basic decimator without the need for
% Matlab's highly optimized signal processing toolbox. This alternate
% decimator's LPF is hamming windowed with cutoff at 0.9pi/M and transition
% band of 0.2pi/M. It also trims the ends of the filtered data (Matlab SPT does
% something clever on the ends to preserve more data.)
%******************************************************************************
function Y=decimate_alt(X,M)

NotRow=~isrow(X); if NotRow, X=X.'; end

%First design a FIR LPF
wc=0.2*pi/M; %cutoff transition band width
N=ceil(3.3/wc); %minimum filter order to meet wc specification
L=N+1; %minimum filter length
if L>size(X,2), L=size(X,2)-1; N=L-1; end %reduce LPF performance if X is too short

%hamming windowed LPF, 0.8*pi/M passband, pi/M stopband 53dB attenuation
h=sinc(0.9*(1/M)*((0:N)-N/2)).*(0.54-0.46*cos(2*pi*(0:N)/N)); h=h/sum(h);

%To decimate apply a LPF and downsample
Y=conv(h,X); %LPF
Y=Y(L:M:end-L); %downsample and trim partially filtered ends

if NotRow, Y=Y.'; end

end
% 
% function h = sinc(t)
% 
% if abs(t) > 1e-12
%     h = sin(pi*t)./t;
% else
%     h = 1;
% end
% 
% end

