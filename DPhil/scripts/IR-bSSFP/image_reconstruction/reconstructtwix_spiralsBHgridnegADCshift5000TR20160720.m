%% Reconstruct Raw Data - Spiral Sampled
clear all

%% 1. Prepare sampled data
%
% set path for code for importing data
addpath(genpath('~/Documents/MATLAB/mapVBVD_20150918'))
%
date = '20160720'; %date of scan
ADCshifts = [0];
ADCshift = 0;
ADCshift_index = find(ADCshifts==(abs(ADCshift)));
scanIndex = ADCshift_index;

MIDstart = 504;
FIDstart = 18946;
% choose data from a specific channel
% simplify by just using the first full sample of k-space
nTRs = 48; %48 is 'fully sampled'


%% 2. Prepare intended k-space location information
%
addpath(genpath('~/Documents/MATLAB/vdspiral'))
addpath(genpath('~/Documents/MATLAB/gridding'))
addpath(genpath('~/Documents/MATLAB/DPhil'))
run('make_spirals')

% k-space coordinates from all the spirals
spiralkspacelocations = [];
for n = 1:nTRs
    tmpspiralkspacelocations = spiralkspace(n,1:nSamplePts);
    spiralkspacelocations = [spiralkspacelocations,tmpspiralkspacelocations];
end
spiralkspacelocations = 0.3*(spiralkspacelocations(:)/(max(spiralkspacelocations(:))));
spiralkspacelocations = spiralkspacelocations';


%% 4. Grid the kspace data
%
spiralkspacelocations_x = real(spiralkspacelocations)';
spiralkspacelocations_y = imag(spiralkspacelocations)';
% determine the density compensation function (DCF)
%   find vertices V and cells C of voronoi diagram
[V,C] = voronoin([spiralkspacelocations_x,spiralkspacelocations_y]);
xVertices = V(:,1);
yVertices = V(:,2);
%   find areas of each cell (DCF)
%for nCell = 1:size(C,1)
%   DCF(nCell) = polyarea(xVertices(C{nCell}),yVertices(C{nCell}));
%end

DCF = [];
for j = 1:length([spiralkspacelocations_x,spiralkspacelocations_y])
    x = V(C{j},1);
    y = V(C{j},2);
    lxy = length(x);
    A = abs(sum(0.5*(x([2:lxy 1])-x(:)).* ...
        (y([2:lxy 1])+y(:))));
    DCF = [DCF A];
end
%removed DCF
DCF(repmat([ones(400,1);zeros(48,1)],48,1)==0)=0; %disregard very high frequency information, to avoid artefacts
DCF = DCF(:)/max(DCF);
% gridding settings
gridsize = 512; %dimension of resulting reconstructed image
kwidth = 1.5; %kernal width
kwidth = 3; %kernal width
overgridfactor = 2;


%%

filename20160720 = (['~/Documents/DPhil/data/TWIX/',date,'/meas_MID',num2str(MIDstart+scanIndex),'_JA_IR_bSSFP_fp_',num2str(abs(ADCshift)),'ADC_FID',num2str(FIDstart+scanIndex)]);
twix_obj20160720 = mapVBVD(filename20160720);
twix_obj20160720.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
% squeeze the data to remove unnecessary dimensions
image_data20160720 = twix_obj20160720.image{''};
TRrange = (4*nTRs)+1:((4*nTRs)+nTRs); %need to be in steady state regime

dat = zeros(32,gridsize,gridsize);
totalReconImage = zeros(1,gridsize*gridsize);
for channel = 1:32
% Read data from file
spiralkspacedata = double(squeeze(image_data20160720(:,channel,TRrange)));

load spiralexampledata
dat(channel,:,:)  = gridkb(spiralkspacelocations,spiralkspacedata,DCF',gridsize,kwidth,overgridfactor);

%% ifft to get images from gridded kspace data
im = fftshift(fft2(fftshift(squeeze(dat(channel,:,:)))));
reconImage(channel,:) = abs(im(:));
end

%combine channels
disp 'Combine Channels...'
for pixel = 1:gridsize*gridsize
    pixel
for channel = 1:32
    totalReconImage(1,pixel) = totalReconImage(1,pixel) + reconImage(channel,pixel)^2;
end
totalReconImage(1,pixel) = sqrt(totalReconImage(1,pixel));
end
disp 'Combine Channels... Complete'
%% 5. Plot Results
figure('Name',['Combined channel data (sum of squares), Grid Recon. ADC shift:',num2str(ADCshift)])
totalReconImage = reshape(totalReconImage,gridsize,gridsize);
imagesc(totalReconImage);
axis square
colorbar


figure('Name',['Grid Recon. ADC shift:',num2str(ADCshift)])
subplot(2,2,1)
plot(spiralkspacelocations(1:448))
title 'spiral k-space locations'
subplot(2,2,2)
plot(spiralkspacedata,'.')
title 'spiral k-space data'
subplot(2,2,3)
plot(squeeze(dat(channel,:,:)),'.')
subplot(2,2,4)
imagesc(reshape(reconImage(channel,:),gridsize,gridsize));
axis square;
title (['Reconstruction'])
