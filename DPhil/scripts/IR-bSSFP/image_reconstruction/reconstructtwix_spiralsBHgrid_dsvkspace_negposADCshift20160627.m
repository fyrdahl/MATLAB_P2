%% Reconstruct Raw Data - Spiral Sampled
clear all

%% 1. Prepare sampled data
%
% set path for code for importing data
addpath(genpath('~/Documents/MATLAB/mapVBVD_20150918'))
%
date = '20160627'; %date of scan

%% 2. Prepare intended k-space location information
%
addpath(genpath('~/Documents/MATLAB/vdspiral'))
addpath(genpath('~/Documents/MATLAB/gridding'))
addpath(genpath('~/Documents/MATLAB/DPhil'))




% read in POET simulation gradients
addpath(genpath('~/Documents/MATLAB/read_dsv'));
% x-axis gradient
dsvGRX_filename = '~/Documents/DPhil/dsv_files/DspData_GRX.dsv';
dsvGRX = Read_dsv(dsvGRX_filename);
%y-gradient
dsvGRY_filename = '~/Documents/DPhil/dsv_files/DspData_GRY.dsv';
dsvGRY = Read_dsv(dsvGRY_filename);
gamma = 4258;
T = 0.00001;

nTotalPts =  567;
nSamplePts = 448;

% convert gradient data to k-space trajectory
GRX = dsv2timecourse(dsvGRX);
GRY = dsv2timecourse(dsvGRY);
spiralkspacelocations_x = cumsum(dsv2timecourse(dsvGRX))*gamma*T;
spiralkspacelocations_y = cumsum(dsv2timecourse(dsvGRY))*gamma*T;
TR = 1000;
%%

nTRs = 48;
%%
spiralkspacelocations_x(1:1839)=[];
spiralkspacelocations_y(1:1839)=[];

%%
for n = 1:nTRs
spiralkspacelocations_x(1+nSamplePts+((n-1)*TR):1000+((n-1)*TR))=0;
spiralkspacelocations_y(1+nSamplePts+((n-1)*TR):1000+((n-1)*TR))=0;

end

%%
spiralkspacelocations_x(spiralkspacelocations_x==0)=[];
spiralkspacelocations_y(spiralkspacelocations_y==0)=[];

spiralkspacelocations_x(1+nTRs*nSamplePts:end) = [];
spiralkspacelocations_y(1+nTRs*nSamplePts:end) = [];

spiralkspacelocations_x = (0.00000001)*spiralkspacelocations_x;
spiralkspacelocations_y = (0.00000001)*spiralkspacelocations_y;

spiralkspacelocations = complex(spiralkspacelocations_x, spiralkspacelocations_y);


ADCshifts = [-100:10:-10,-7,-3,0,3,7,10:10:100];
ADCshift = 0% [-10, -7, -3, 0, 3, 7, 10] %ADCshifts;
ADCshift_index = find(ADCshifts==ADCshift);
scanIndex = ADCshift_index;


channel =1
% read data from file
filename20160627 = (['~/Documents/DPhil/data/TWIX/',date,'/meas_MID',num2str(10+scanIndex),'_JA_IR_bSSFP_fp_NOINV_',num2str(ADCshift),'_FID',num2str(16618+scanIndex)]);
twix_obj20160627 = mapVBVD(filename20160627);
twix_obj20160627.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
% squeeze the data to remove unnecessary dimensions
image_data20160627 = twix_obj20160627.image{''};
%image_data(440:end,:,:) = 0;
spiralkspacedata = double(squeeze(image_data20160627(:,channel,((3*nTRs)+1):((3*nTRs)+nTRs) )));
figure
for n = 1:48
plot(spiralkspacedata)
end
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
DCF(repmat([ones(440,1);zeros(8,1)],48,1)==0)=DCF(439);
DCF = DCF(:)/max(DCF);
% gridding settings
gridsize = 256; %dimension of resulting reconstructed image
kwidth = 1.5; %kernal width
overgridfactor = 2;
load spiralexampledata
[dat]  = gridkb(spiralkspacelocations,spiralkspacedata,DCF',gridsize,kwidth,overgridfactor);
%% 4. Plot Results
%% 4. Plot Results
figure('Name',[date,' Grid Recon. ADC shift:',num2str(ADCshift),' channel',num2str(channel)])
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


