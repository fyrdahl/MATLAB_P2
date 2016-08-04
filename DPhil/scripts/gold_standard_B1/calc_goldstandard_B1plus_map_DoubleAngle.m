%% Mapping the B1+ Profile - Double-Angle Method
% see:
% http://web.stanford.edu/class/rad229/Notes/4d-Mapping.pdf
% Neeb, MRM 2008
% Morrell, Physics in medicine and biology 2010

workingdir = '~/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze

savingdir = '~/Documents/DPhil';
addpath(genpath(savingdir));

% Load images
image1filepath = ['~/Documents/DPhil/data/NIfTI/20160503_MR_PHYSICS_462/1_002_gre_TR3000_FA45_20160503/20160503_002_gre_TR3000_FA45'];
B1plusimage1 = read_avw(image1filepath);
image2filepath = ['~/Documents/DPhil/data/NIfTI/20160503_MR_PHYSICS_462/1_003_gre_TR3000_FA90_20160503/20160503_003_gre_TR3000_FA90'];
B1plusimage2 = read_avw(image2filepath);


% effective flip angle (see: Neeb 2008 and Hargreaves pdf slides)
effAlpha = acos(B1plusimage2./(2*B1plusimage1));


figure, 
imagesc(real(effAlpha))
colorbar
caxis([0.5 1.2])
title 'B1+ Efficiency - Double Angle method'