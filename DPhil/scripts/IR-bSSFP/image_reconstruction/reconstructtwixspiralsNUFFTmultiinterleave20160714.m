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
run('make_spirals')

% k-space coordinates from all the spirals
spiralkspacelocations = [];
for n = 1:nTRs
    tmpspiralkspacelocations = spiralkspace(n,1:nSamplePts);
    spiralkspacelocations = [spiralkspacelocations,tmpspiralkspacelocations];
end
% k-space location radii should not exceed 0.5
spiralkspacelocations = 0.3*(spiralkspacelocations(:)/(max(spiralkspacelocations(:))));
spiralkspacelocations = spiralkspacelocations';
%% Load raw data
disp('Load Raw Data...')
ADCshifts = 0:10:110;

date = '20160714'; %date of scan
ADCshifts = [2:2:10];
ADCshift = 4;
ADCshift_index = find(ADCshifts==ADCshift);
scanIndex = ADCshift_index;

% choose data from a specific channel
% simplify by just using the first full sample of k-space
nTRs = 48; %48 is 'fully sampled'
channel = 1;
MIDstart = 576;
FIDstart = 17955;
% Read data from file
filename20160714 = (['~/Documents/DPhil/data/TWIX/',date,'/meas_MID',num2str(MIDstart+scanIndex),'_JA_IR_bSSFP_fp_neg',num2str(ADCshift),'ADC_FID',num2str(FIDstart+scanIndex)]);

twix_obj20160714 = mapVBVD(filename20160714);
twix_obj20160714.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
% squeeze the data to remove unnecessary dimensions
image_data20160714 = twix_obj20160714.image{''};
disp('Load Raw Data...Finished')


%% NUFFT
%
disp('Prepare Data Vectors for NUFFT...')
nTRs = 48;
nSamplePts = 488;
leaf_values = [];
spiralkspacelocations
for n = (1+2*nTRs):nTRs+2*nTRs
 tmpleaf_values = cast(squeeze(image_data20160714(:,channel,n)),'double') %measured k-space data values
 leaf_values = [leaf_values;tmpleaf_values];
end
spiralkspacelocations = complex(real(spiralkspacelocations)/max(real(spiralkspacelocations)),imag(spiralkspacelocations)/max(imag(spiralkspacelocations))) ;
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
