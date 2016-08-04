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

%% Load raw data
disp('Load Raw Data...')
ADCshifts = 0:10:110;

date = '20160416';
nChannel = 1;
ADCshift = ADCshifts(1);
filename = (['/Users/jallen/Documents/DPhil/data/TWIX/',date,'/meas_MID7',num2str(28+(ADCshifts(1+(ADCshift*0.1))*0.1)),'_JA_IR_bSSFP_fp_tr10_fa10_noinv_shift',num2str(ADCshifts(1+(ADCshift*0.1))),'_FID',num2str(9497+(ADCshifts(1+(ADCshift*0.1))*0.1))]);
twix_obj = mapVBVD(filename);
twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
image_data = twix_obj.image{''};
disp('Load Raw Data...Finished')

%% NUFFT
%
disp('Prepare Data Vectors for NUFFT...')
numInterleaves = 48;
for n = 1:numInterleaves
 leaf_values((n-1)*nSamplePts+1:n*nSamplePts) = cast(squeeze(image_data(:,nChannel,n+(2*48))),'double'); %measured k-space data values
 samplelocations(1,(n-1)*nSamplePts+1:n*nSamplePts) = spiralkspace(n,1:nSamplePts); %desired k-space sample positions
end
% magnitude of radii must be < 5
samplelocations = complex(0.5*real(samplelocations)/max(real(samplelocations)),0.5*imag(samplelocations)/max(imag(samplelocations))) ;
radii = sqrt((real(samplelocations)).^2 + (imag(samplelocations)).^2);
%
MLfile = [real(samplelocations); imag(samplelocations)]';
MLvec_of_samples = leaf_values;
fftvalues=[MLfile MLvec_of_samples'];
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
