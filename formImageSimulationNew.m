    function data= formImageSimulationNew(elev,Nx,Ny,azimuthVals,...
        phaseHistory,f,Lsize_x,Lsize_y,Nfft,taper_flag,AntX,AntY,AntZ,R0,x_vec,y_vec,x0,y0)

% Define image parameters here
data.Wx = Lsize_x;           % Scene extent x (m)
data.Wy = Lsize_y;           % Scene extent y (m)
data.Nfft = Nfft;   % Number of samples in FFT
data.R0 =  R0;
data.Nx = Nx;          % Number of samples in x direction
data.Ny = Ny;          % Number of samples in y direction
data.x0 = x0;            % Center of image scene in x direction (m)
data.y0 = y0;            % Center of image scene in y direction (m)
dyn_range = 70;         % dB of dynamic range to display

% INPUT PARAMETERS END HERE %

% Update the phase history

data.phdata = phaseHistory;
% Update other parameters needed for imaging
data.AntAzim = azimuthVals;
data.AntAz = azimuthVals*pi/180; % radians
data.AntElev =elev;
data.freq = f ;

% Calculate the minimum frequency for each pulse (Hz)
data.minF = min(data.freq)*ones(size(data.AntAzim));

% Calculate the frequency step size (Hz)
data.deltaF = diff(data.freq(1:2));

% Determine the number of pulses and the samples per pulse
[data.K,data.Np] = size(data.phdata);
data.AntX = AntX;
data.AntY = AntY;
data.AntZ = AntZ;

% Add a hamming taper to the data if desired
if taper_flag
    data.phdata = data.phdata .* (hamming(data.K)*hamming(data.Np)');
  %  data.phdata = data.phdata .* (taylorwin(data.K,4,-35)*taylorwin(data.Np,4,-35)');
end

% Setup imaging grid
data.x_vec = x_vec;
data.y_vec = y_vec;

% data.x_vec = linspace(data.x0 - data.Wx/2, data.x0 + data.Wx/2, data.Nx);
% data.y_vec = linspace(data.y0 - data.Wy/2, data.y0 + data.Wy/2, data.Ny);
[data.x_mat,data.y_mat] = meshgrid(data.x_vec,data.y_vec);
data.z_mat = zeros(size(data.x_mat));

% Call the backprojection function with the appropriate inputs
data = bpBasic_mod(data);

% Display the image
% figure
% imagesc(data.x_vec,data.y_vec,20*log10(abs((data.im_final))./...
%     max(max(abs(data.im_final)))),[-dyn_range 0])
% colormap jet
% axis xy image;
% set(gca,'XTick',-5:5,'YTick',-5:5);
% h = xlabel('x (m)');
% set(h,'FontSize',14,'FontWeight','Bold');
% h = ylabel('y (m)');
% set(h,'FontSize',14,'FontWeight','Bold');
% colorbar
% set(gca,'FontSize',14,'FontWeight','Bold');
% print -deps2 /ssip2/lgorham/SPIE10/fig/CVdomesBPA.eps