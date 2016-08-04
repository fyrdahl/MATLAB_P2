%% Reconstruct Raw Data - Spiral Sampled
clear all

%% 1. Prepare sampled data
%
% set path for code for importing data
addpath(genpath('~/Documents/MATLAB/mapVBVD_20150918'))
%
date = '20160415'; %date of scan

%% 2. Prepare intended k-space location information
%
addpath(genpath('~/Documents/MATLAB/vdspiral'))
addpath(genpath('~/Documents/MATLAB/gridding'))
addpath(genpath('~/Documents/MATLAB/DPhil'))
run('make_spirals')
% k-space coordinates from all the spirals
spiralkspacelocations = [];

channel = 1;
% choose data from a specific channel
% simplify by just using the first full sample of k-space
nTRs = 48; %48 is 'fully sampled'

for n = 1:nTRs
tmpspiralkspacelocations = spiralkspace(n,1:nSamplePts);
spiralkspacelocations = [spiralkspacelocations,tmpspiralkspacelocations];
end
% k-space location radii should not exceed 0.5
spiralkspacelocations = 0.3*(spiralkspacelocations(:)/(max(spiralkspacelocations(:))));
spiralkspacelocations = spiralkspacelocations';


ADCshifts = 0:10:110;
for ADCshift_index = 1%:numel(ADCshifts);
ADCshift = ADCshifts(ADCshift_index);
scanIndex = ADCshift_index;

% read data from file
filename = (['~/Documents/DPhil/data/TWIX/',date,'/meas_MID',num2str(727+scanIndex),'_JA_IR_bSSFP_fp_tr10_fa10_noinv_shift',num2str(ADCshift),'_FID',num2str(9496+scanIndex)]);
twix_obj = mapVBVD(filename);
twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory 
% squeeze the data to remove unnecessary dimensions
image_data = twix_obj.image{''};
image_data(440:end,:,:) = 0;

spiralkspacedata = double(squeeze(image_data(:,channel,1:nTRs)));


%% 3. Grid the kspace data
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
DCF(repmat([ones(440,1);zeros(8,1)],48,1)==0)=0;
DCF = DCF(:)/max(DCF);
% gridding settings
gridsize = 512; %dimension of resulting reconstructed image
kwidth = 1.5; %kernal width
overgridfactor = 2;
load spiralexampledata
[dat]  = gridkb(spiralkspacelocations,spiralkspacedata,DCF',gridsize,kwidth,overgridfactor);
%% 4. Plot Results
figure('Name',['Grid Recon. ADC shift:',num2str(ADCshift)])
subplot(2,2,1)
plot(spiralkspacelocations)
title 'spiral k-space locations'
subplot(2,2,2)
plot(spiralkspacedata,'.')
title 'spiral k-space data'
subplot(2,2,3)
plot(dat,'.')
subplot(2,2,4)
im = fftshift(fft2(fftshift(dat)));
im = abs(im);
imagesc(im);
axis square;
title (['Reconstruction'])
end