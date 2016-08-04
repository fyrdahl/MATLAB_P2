%% Reconstruct Raw Data - Spiral Sampled
clear all

%% 1. Prepare sampled data
%
% set path for code for importing data
addpath(genpath('~/Documents/MATLAB/mapVBVD_20150918'))
%
date = '20160629'; %date of scan

%% 2. Prepare intended k-space location information
%
addpath(genpath('~/Documents/MATLAB/vdspiral'))
addpath(genpath('~/Documents/MATLAB/gridding'))
addpath(genpath('~/Documents/MATLAB/DPhil'))

% read in POET simulation gradients
addpath(genpath('~/Documents/MATLAB/read_dsv'));
% x-axis gradient
dsvGRX_filename = '~/Documents/DPhil/dsv_files/20160629/DspData_GRX.dsv';
dsvGRX = Read_dsv(dsvGRX_filename);
%y-gradient
dsvGRY_filename = '~/Documents/DPhil/dsv_files/20160629/DspData_GRY.dsv';
dsvGRY = Read_dsv(dsvGRY_filename);
%
gamma = 4258;
T = 0.00001;
nTRs = 48;
nTotalPts =  594;
nSamplePts = 488;
TR = 1000;

% convert gradient data to k-space trajectory
GRX = dsv2timecourse(dsvGRX);
GRY = dsv2timecourse(dsvGRY);
spiralkspacelocations_x = cumsum(dsv2timecourse(dsvGRX))*gamma*T;
spiralkspacelocations_y = cumsum(dsv2timecourse(dsvGRY))*gamma*T;
%
% Remove dead time between readouts
initialdeadtime = 540;
spiralkspacelocations_x(1:initialdeadtime)=[];
spiralkspacelocations_y(1:initialdeadtime)=[];
for n = 1:nTRs
    spiralkspacelocations_x(1+nSamplePts+((n-1)*TR):1000+((n-1)*TR))=0;
    spiralkspacelocations_y(1+nSamplePts+((n-1)*TR):1000+((n-1)*TR))=0;
end
spiralkspacelocations_x(spiralkspacelocations_x==0)=[];
spiralkspacelocations_y(spiralkspacelocations_y==0)=[];
spiralkspacelocations_x(1+nTRs*nSamplePts:end) = [];
spiralkspacelocations_y(1+nTRs*nSamplePts:end) = [];
% Re-scale
spiralkspacelocations_x = (0.00000001)*spiralkspacelocations_x;
spiralkspacelocations_y = (0.00000001)*spiralkspacelocations_y;
%
spiralkspacelocations = complex(spiralkspacelocations_x, spiralkspacelocations_y);


%% Read data
channel = 1;
ADCshifts = [-100:10:-10,-7,-4:1:4,7,10:10:100];
%for ADCshift = [-7, -4:1:4, 7] %0 ;%ADCshifts;
ADCshift = 0 
ADCshift_index = find(ADCshifts==ADCshift);
scanIndex = ADCshift_index;
%
% read data from file
filename20160629 = (['~/Documents/DPhil/data/TWIX/',date,'/meas_MID',num2str(269+scanIndex),'_JA_IR_bSSFP_fp_NOINV_',num2str(ADCshift),'_FID',num2str(16876+scanIndex)]);
twix_obj20160629 = mapVBVD(filename20160629);
twix_obj20160629.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
% squeeze the data to remove unnecessary dimensions
image_data20160629 = twix_obj20160629.image{''};
%image_data(470:end,:,:) = 0;
spiralkspacedata = double(squeeze(image_data20160629(:,channel,(1:48)+(2*48)) ));
figure
for n = 1:48
plot(spiralkspacedata(:,n))
hold on
title([date,' ADCshift:',num2str(ADCshift)])
end
%end

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
DCF(repmat([ones(470,1);zeros(18,1)],48,1)==0)=DCF(469); % CHANGE 440 TO RELEVANT INDEX FOR THIS DATA
DCF = DCF(:)/max(DCF);
% gridding settings
gridsize = 256; %dimension of resulting reconstructed image
kwidth = 1.5; %kernal width
overgridfactor = 2;
load spiralexampledata
[dat]  = gridkb(spiralkspacelocations,spiralkspacedata,DCF',gridsize,kwidth,overgridfactor);
%% 4. Plot Results
figure('Name',[date,' Grid Recon. ADC shift:',num2str(ADCshift)])
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
