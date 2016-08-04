%% Reconstruct Raw Data using NUFFT - Spiral Sampled
%
%
%% 1. File and Path Preparation
addpath('/Users/jallen/Documents/MATLAB/irt',...
    '/Users/jallen/Documents/MATLAB/vdspiral')
addpath(genpath('/Users/jallen/Documents/MATLAB/DPhil'),...
    genpath('/Users/jallen/Documents/MATLAB/mapVBVD_20150918'))
run('/Users/jallen/Documents/MATLAB/irt/setup')

%% 2. Load raw data
disp('Load Raw Data...')
%
ADCshifts = 0:10:110;
nChannel = 1;
%
for ADCshift = ADCshifts(1)
    filename = (['/Users/jallen/Documents/DPhil/data/TWIX/20160416/meas_MID7',num2str(28+(ADCshifts(1+(ADCshift*0.1))*0.1)),'_JA_IR_bSSFP_fp_tr10_fa10_noinv_shift',num2str(ADCshifts(1+(ADCshift*0.1))),'_FID',num2str(9497+(ADCshifts(1+(ADCshift*0.1))*0.1))]);
    twix_obj = mapVBVD(filename);
    twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
    image_data = twix_obj.image{''};
    kimage = squeeze(image_data(:,nChannel,:));
end
disp('Load Raw Data...Finished')
%
%% 3. Spiral Reconstruction
disp('Prepare Spiral Trajectories...')
run('make_spirals')
nADCshift = 1+(ADCshift*0.1);
% k-space coordinates from all the spirals
allKTrajs = zeros(nSamplePts,2);
tmpksamples = zeros(nSamplePts,2);
for nADCshift = 1:nInterleaves
    allKTrajs((1 + ((nADCshift-1)*nSamplePts)):(nSamplePts + ((nADCshift-1)*nSamplePts)),1) = cumsum(spiralgradients(nADCshift,1:nSamplePts))*gamma*T;
    allKTrajs((1 + ((nADCshift-1)*nSamplePts)):(nSamplePts + ((nADCshift-1)*nSamplePts)),2) = cumsum(spiralgradients(nADCshift+nInterleaves,1:nSamplePts))*gamma*T;
end
disp('Prepare Spiral Trajectories...Finished')
%% 4. NUFFT
disp('Prepare Data Vectors for NUFFT...')
nleaf = 1;
single_leaf_values=kimage(:,nleaf);
firstspiralkspace_locations=spiralkspace(nleaf,1:448);
vec_of_samples = 0.5*firstspiralkspace_locations/max(max(firstspiralkspace_locations));
fft_values=[real(single_leaf_values) imag(single_leaf_values) vec_of_samples'];
disp('Prepare Data Vectors for NUFFT...Finished')
disp('Perform Reconstruction...')
[ reconstructed_image ] = non_uniform_fessler_ifft2c(fft_values,256); % Lior's Function, using Fessler's (Michigan) Code
disp('Perform Reconstruction...Finished')
figure, imagesc(abs(reconstructed_image)), title('Single Interleave: Reconstructed Image')
dt = datetime('now');
%save([savingdir,'/MAT-files/images/reconstructed_image',datestr(dt),'.mat'],'reconstructed_image')
%matlab2tikz([savingdir,'/MAT-files/images/reconstructed_image',datestr(dt),'.tex'])
