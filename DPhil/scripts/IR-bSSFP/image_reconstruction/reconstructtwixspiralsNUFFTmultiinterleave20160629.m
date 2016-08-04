%% Reconstruct Raw Data using NUFFT - Spiral Sampled
%
%
%% 1. Set up paths
%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/DPhil';
workingdir = '/Users/jallen/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze
%
savingdir = '/Users/jallen/Documents/DPhil';
addpath(genpath(savingdir));

% If working on jalapeno00, uncomment the following lines:
% addpath(genpath('/Applications/fsl/'))
% addpath(genpath('/usr/local/fsl/bin'))
% addpath(genpath('/opt/fmrib/fsl/etc/matlab'))
addpath('/Users/jallen/Documents/MATLAB/irt',...
    '/Users/jallen/Documents/MATLAB/vdspiral')
addpath(genpath('/Users/jallen/Documents/MATLAB/DPhil'),...
    genpath('/Users/jallen/Documents/MATLAB/mapVBVD_20150918'))
run('/Users/jallen/Documents/MATLAB/irt/setup')

%% Spiral Reconstruction
disp('Prepare Spiral Trajectories...')
run('make_spirals')
disp('Prepare Spiral Trajectories...Finished')
%% Spiral Locations
% convert gradient data to k-space trajectory
%
addpath(genpath('~/Documents/MATLAB/vdspiral'))
addpath(genpath('~/Documents/MATLAB/gridding'))
addpath(genpath('~/Documents/MATLAB/DPhil'))
%
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

%% Load raw data
disp('Load Raw Data...')
ADCshifts = 0:10:110;

date = '20160629';
channel = 1;
ADCshifts = [-100:10:-10,-7,-4:1:4,7,10:10:100];
%for ADCshift = [-7, -4:1:4, 7] %0 ;%ADCshifts;
ADCshift = 0 
ADCshift_index = find(ADCshifts==ADCshift);
scanIndex = ADCshift_index;
%
filename = (['~/Documents/DPhil/data/TWIX/',date,'/meas_MID',num2str(269+scanIndex),'_JA_IR_bSSFP_fp_NOINV_',num2str(ADCshift),'_FID',num2str(16876+scanIndex)]);
twix_obj = mapVBVD(filename);
twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
image_data = twix_obj.image{''};
disp('Load Raw Data...Finished')


%% NUFFT
%
disp('Prepare Data Vectors for NUFFT...')
nTRs = 48;
nSamplePts = 488;
leaf_values = [];
spiralkspacelocations
for n = (1+2*nTRs):nTRs+2*nTRs
 tmpleaf_values = cast(squeeze(image_data(:,channel,n)),'double') %measured k-space data values
 leaf_values = [leaf_values;tmpleaf_values];
end
spiralkspacelocations = complex(0.5*real(spiralkspacelocations)/max(real(spiralkspacelocations)),0.5*imag(spiralkspacelocations)/max(imag(spiralkspacelocations))) ;
%%
MLfile = [real(spiralkspacelocations); imag(spiralkspacelocations)]';
MLvec_of_samples = leaf_values;
fftvalues=[MLfile MLvec_of_samples];
disp('Prepare Data Vectors for NUFFT...Finished')
%

disp('Perform NUFFT...')
[ reconstructed_image ] = non_uniform_fessler_ifft2c(fftvalues,256); % Lior's Function, using Fessler's (Michigan) Code
disp('Perform NUFFT...Finished')
dt = datetime('now');
abs_reconstructed_image = abs(reconstructed_image);
figure, imagesc(abs_reconstructed_image), axis square
%% save image
%save([savingdir,'/MAT-files/images/reconstructed_imageMI',datestr(dt),'.mat'],'abs_reconstructed_image')
%matlab2tikz([savingdir,'/MAT-files/images/reconstructed_imageMI',datestr(dt),'.tex'])
