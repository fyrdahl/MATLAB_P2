%% Reconstruct Raw Data - Spiral Sampled
clear all
addpath(genpath('/Users/jallen/Documents/MATLAB/mapVBVD_20150918'))

ADCshifts = [0:10:110];
ADCshift = ADCshifts(1);
dt =1;

% 1. Prepare sampled data
% choose data from a specific channel
% simplify by just using the first full sample of k-space
channel = 1;
nTRs = 48; %48 is 'fully sampled'
% read data from file
filename = (['/Users/jallen/Documents/DPhil/data/TWIX/20160415/meas_MID7',num2str(28+(ADCshift/10)),...
    '_JA_IR_bSSFP_fp_tr10_fa10_noinv_shift',num2str(ADCshift),'_FID',num2str(9497+(ADCshift/10))])
twix_obj = mapVBVD(filename);
twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
% squeeze the data to remove unnecessary dimensions
image_data = twix_obj.image{''};

spiralkspacedata = double(squeeze(image_data(:,channel,1:nTRs)));
for n=1:nTRs
    spiralkspacedata(:,n) = fftshift(squeeze(spiralkspacedata(:,n)),2);
end
% 2. Prepare intended k-space locations
%
addpath(genpath('~/Documents/MATLAB/vdspiral'))
addpath(genpath('~/Documents/MATLAB/gridding'))
addpath(genpath('~/Documents/MATLAB/DPhil'))
run('make_spirals')

spiralkspacelocations = [];
for n = 1:nTRs
    % optional truncation of Spiral Locations
    if ADCshift < 0; %start ADC before gradients
        tmpspiralkspacelocations = [zeros(1,abs(ADCshift*dt)),spiralkspace(n,1:(nSamplePts-abs(ADCshift*dt)) )];
    end
    if ADCshift > 0; %start ADC after gradients
        tmpspiralkspacelocations = spiralkspace(n,(1+abs(ADCshift*dt)):(1+abs(ADCshift*dt)+nSamplePts-1) );
    end
    if ADCshift == 0; % start ADC exactly same time as readout gradients
        tmpspiralkspacelocations = spiralkspace(n,1:nSamplePts);
    end
    spiralkspacelocations = [spiralkspacelocations,tmpspiralkspacelocations];
    
end

%% k-space location radii should not exceed 0.5
spiralkspacelocations = pi*(1/10)*spiralkspacelocations/(max(spiralkspacelocations));
%spiralkspacelocations = spiralkspacelocations(:)';




%% 3. Grid the kspace data
%
spiralkspacelocations_x = real(spiralkspacelocations)';
spiralkspacelocations_y = imag(spiralkspacelocations)';

% determine the density compensation function (DCF)
%   find vertices V and cells C of voronoi diagram
[V,C] = voronoin([spiralkspacelocations_x,spiralkspacelocations_y]);
xVertices = V(:,1);
yVertices = V(:,2);
%   find areas of each cell
for nCell = 1:size(C,1)
    DCF(nCell) = polyarea(xVertices(C{nCell}),yVertices(C{nCell}));
end
%   remove DCF values for outer k-space
DCF(repmat([ones(440,1);zeros(8,1)],48,1)==0)=0;
DCF = DCF(:)/max(DCF);

% gridding settings
gridsize = 256; %dimension of resulting reconstructed image
kwidth = 1.5; %kernal width
overgridfactor = 2;
load spiralexampledata

% grid data
[dat]  = gridkb(spiralkspacelocations,spiralkspacedata,DCF',gridsize,kwidth,overgridfactor);

%% 4. Plot Results
figure,
subplot(2,2,1)
plot(spiralkspacelocations)
title 'spiral k-space locations'
subplot(2,2,2)
plot(spiralkspacedata,'.')
title 'spiral k-space data'
subplot(2,2,3)
plot(dat,'.')
subplot(2,2,4)
im = ifftshift(fft2(((dat))));
im = abs(im);
imagesc(im);
axis square;

