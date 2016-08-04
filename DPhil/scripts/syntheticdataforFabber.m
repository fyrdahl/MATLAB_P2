%% ANALYSIS PROJECT - Testing FABBER
% 7th June 2016

%% Notes:
% trying terminal line options:
%./fabber --output=data_out --model=VFA --TR=500 --FAvals=fabberFAvals.txt
%--method=vb --noise=white --mask=/Users/jallen/Documents/DPhil/data/NIfTI/20150714_MR_PHYSICS_338/TI_group/20150714_002_IR_ep2d_se_TI_35ms.nii
%--data1=/Users/jallen/Documents/DPhil/data/NIfTI/20150714_MR_PHYSICS_338/TI_group/20150714_002_IR_ep2d_se_TI_35ms.nii
%--data2=/Users/jallen/Documents/DPhil/data/NIfTI/20150714_MR_PHYSICS_338/TI_group/20150714_003_IR_ep2d_se_TI_85ms.nii


%% generate synthetic data
FA1 = 45;
FA2 = 90;
TE = 500;

% synthetic T1 map
T1map = phantom*1000 + 1000;

% post-FA
image1 = T1map*(FA1/90);
image2 = T1map*(FA2/90);

%exponential decay
image1d = image1.*exp(-TE./T1map);
image2d = image2.*exp(-TE./T1map);

image1nii=make_nii(image1d)
save_nii(image1nii,'/Users/jallen/Documents/fabber/image1nii.nii')

image2nii=make_nii(image2d)
save_nii(image2nii,'/Users/jallen/Documents/fabber/image2nii.nii')


figure, 
subplot 221
imagesc(image1)
caxis([0 2000])
subplot 222
imagesc(image2)
caxis([0 2000])
subplot 223
imagesc(image1d)
caxis([0 2000])
subplot 224
imagesc(image2d)
caxis([0 2000])