%% Reconstruct Cartesian Sampled Localiser Data
clear all
%
filename = '/Users/jallen/Documents/DPhil/data/TWIX/20160415/meas_MID727_localizer_FID9496';

twix_obj = mapVBVD(filename);
image_data = twix_obj.image{''};

totalImage = zeros(size(squeeze(image_data(:,1,:,1,1))));

for n = 1:32
kimage = squeeze(image_data(:,n,:,1,1));
ifftkimage(:,:,n) = ifftshift(ifft2(kimage));
totalImage = totalImage + squeeze(ifftkimage(:,:,n));
end

%%
nCh = 1;
figure
subplot(1,2,1)
imagesc(abs(squeeze(ifftkimage(:,:,nCh))))
title(['Channel ',num2str(nCh)])
subplot(1,2,2)
imagesc(abs(totalImage))
title 'Addition of all channels'


%% Spiral Reconstruction

%% NUFFT initiation variables
%|	om [M d]	"digital" frequencies in radians
% k-space coordinates from all the spirals
allKSamples = zeros(nSamplePts,2);
tmpksamples = zeros(nSamplePts,2);
for n = 1:nInterleaves
   allKSamples((1 + ((n-1)*nSamplePts)):(nSamplePts + ((n-1)*nSamplePts)),1) = cumsum(spiralgradients(n,1:nSamplePts))*gamma*T;
   allKSamples((1 + ((n-1)*nSamplePts)):(nSamplePts + ((n-1)*nSamplePts)),2) = cumsum(spiralgradients(n+nInterleaves,1:nSamplePts))*gamma*T;
end
om = allKSamples;
%|	Nd [d]		image dimensions (N1,N2,...,Nd)
Nd = size(kimage); % 
%|	Jd [d]		# of neighbors used (in each direction)
Jd = [6 6]; % same used in nufft_example.m from Jeff Fessler
%|	Kd [d]		FFT sizes (should be >= N1,N2,...)
Kd = 2*Nd; % A factor of 2 was used in nufft_example.m

%%
st = nufft_init(allKSamples, size(kimage), Jd, Kd);
testImage = nufft(kimage,st);