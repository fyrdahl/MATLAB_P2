%% Reconstruct Raw Data using NUFFT - Spiral Sampled

%% File and Path Preparation
addpath('/Users/jallen/Documents/MATLAB/irt',...
    '/Users/jallen/Documents/MATLAB/vdspiral')
addpath(genpath('/Users/jallen/Documents/MATLAB/DPhil'))
run('/Users/jallen/Documents/MATLAB/irt/setup')

%% Load raw data
disp('Load Raw Data...')
ADCshifts = [0:10:110];
nChannel = 1;
for ADCshift = ADCshifts(5)
filename = (['/Users/jallen/Documents/DPhil/data/TWIX/20160415/meas_MID7',num2str(28+(ADCshifts(1+(ADCshift*0.1))*0.1)),'_JA_IR_bSSFP_fp_tr10_fa10_noinv_shift',num2str(ADCshifts(1+(ADCshift*0.1))),'_FID',num2str(9497+(ADCshifts(1+(ADCshift*0.1))*0.1))]);
twix_obj = mapVBVD(filename);
twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory 
image_data = twix_obj.image{''};
kimage = squeeze(image_data(:,nChannel,:));
end
disp('Load Raw Data...Finished')
%%
%{
% Plot k-space samples for all ADC shifts
figure
for nADCshift = 1:numel(ADCshifts)
subplot(3,4,nADCshift)
plot(squeeze(kspace(:,:,nADCshift,nChannel)),'*')
drawnow;
end
%}
%% Spiral Reconstruction
disp('Prepare Spiral Trajectories...')
run('make_spirals')
nADCshift = 1+(ADCshift*0.1);
% NUFFT initiation variables
%|	om [M d]	"digital" frequencies in radians
% k-space coordinates from all the spirals
allKTrajs = zeros(nSamplePts,2);
tmpksamples = zeros(nSamplePts,2);
for nADCshift = 1:nInterleaves
   allKTrajs((1 + ((nADCshift-1)*nSamplePts)):(nSamplePts + ((nADCshift-1)*nSamplePts)),1) = cumsum(spiralgradients(nADCshift,1:nSamplePts))*gamma*T;
   allKTrajs((1 + ((nADCshift-1)*nSamplePts)):(nSamplePts + ((nADCshift-1)*nSamplePts)),2) = cumsum(spiralgradients(nADCshift+nInterleaves,1:nSamplePts))*gamma*T;
end
disp('Prepare Spiral Trajectories...Finished')
%% NUFFT
kx = cast(real(kimage(:,1)),'double');
ky = cast(imag(kimage(:,1)),'double');
om = allKTrajs;
%
%|	Nd [d]		image dimensions (N1,N2,...,Nd)
Nd = [256 256]; % 
%|	Jd [d]		# of neighbors used (in each direction)
Jd = [6 6]; % same used in nufft_example.m from Jeff Fessler
%|	Kd [d]		FFT sizes (should be >= N1,N2,...)
Kd = 2*Nd; % A factor of 2 was used in nufft_example.m
%
% Initialise structure for NUFFT
st = nufft_init(allKTrajs, size(kimage), Jd, Kd);
% nufft returns a spectra 'X'
X = nufft(kimage,st);
%
[oo1 oo2] = ndgrid(	2*pi*([0:Nd(1)-1]/Nd(1) - 0.5), ...
2*pi*([0:Nd(2)-1]/Nd(2) - 0.5));
% omega is taken from the trajectory. A 2xnPts matrix varying between -pi and pi
omega = 2*pi*[allKTrajs(:,1) allKTrajs(:,2)] / Nd(1);
omega = cast(omega,'double');
% yi is complex vector of data
yi(:,1) = real(reshape(kimage(:,1:48),1,size(omega,1)));
yi(:,2) = imag(reshape(kimage(:,1:48),1,size(omega,1)));
yi = cast(yi,'double');
%
yd_g = griddata(omega(:,1), omega(:,2), yi(1:size(omega,1),1), oo1, oo2, 'cubic');
yd_g(isnan(yd_g)) = 0;

wi = abs(omega(:,1) + 1i * omega(:,2));
xcp = X' * (wi .* yi);

%%
xg = ifft2(fftshift(yd_g));
figure, imagesc(abs(xg))

Xg = ifft2(fftshift(X));
figure, imagesc(abs(Xg))

