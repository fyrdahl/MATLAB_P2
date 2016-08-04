%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/DPhil';
workingdir = '/Users/jallen/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze

savingdir = '/Users/jallen/Documents/DPhil';
addpath(genpath(savingdir));

% If working on jalapeno00, uncomment the following lines:
% addpath(genpath('/Applications/fsl/'))
% addpath(genpath('/usr/local/fsl/bin'))
% addpath(genpath('/opt/fmrib/fsl/etc/matlab'))

%% 2. Load image

filepath = [savingdir,'/data/NIfTI/20160426_MR_PHYSICS_459/1_017_gre_field_mapping_3000TR_50TE_20160426/20160426_017_gre_field_mapping_3000TR_50TE']
FPimages = read_avw(filepath);

FMimage1 = squeeze(FPimages(:,:,9));
FMimage2 = squeeze(FPimages(:,:,10));
FMimage3 = squeeze(FPimages(:,:,11));
FMimage4 = squeeze(FPimages(:,:,12));

%% 3. Image Ratios

ratio = FMimage1./FMimage2;
1./(2*ratio)
alpha = acos(1./(2*ratio));


%% Plot Results


figure, 
subplot(1,2,1)
imagesc(ratio)
subplot(1,2,2)
imagesc(real(alpha))



%% original images
figure, 
subplot(2,2,1)
imagesc(FMimage1)
title 'slice 9'
subplot(2,2,2)
imagesc(FMimage2)
title 'slice 10'
subplot(2,2,3)
imagesc(FMimage3)
title 'slice 11'
subplot(2,2,4)
imagesc(FMimage4)
title 'slice 12'