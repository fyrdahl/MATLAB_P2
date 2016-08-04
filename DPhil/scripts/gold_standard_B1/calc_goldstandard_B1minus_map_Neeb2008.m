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
image1filepath = ['~/Documents/DPhil/data/NIfTI/20160503_MR_PHYSICS_462/1_004_gre_TR3000_FA90_BODYREC_20160503/20160503_004_gre_TR3000_FA90_BODYREC'];
B1minusimageBR = read_avw(image1filepath);
image2filepath = ['~/Documents/DPhil/data/NIfTI/20160503_MR_PHYSICS_462/1_005_gre_TR3000_FA90_20160503/20160503_005_gre_TR3000_FA90'];
B1minusimageHR = read_avw(image2filepath);


%% Calculate B1- Map
% effective flip angle (see: Neeb 2008 and Hargreaves pdf slides)

%Sensitivity = image1/(image2*(effectiveFA/nominalFA))
B1minusprofileNeeb = B1minusimageBR./(B1minusimageHR.*((real(effAlpha).*90)./90)); %Neeb2008



%% Plot Results
figure,
subplot(2,3,1)
imagesc(B1minusimageBR)
title 'GRE TR3000 FA90 Body Receive'
subplot(2,3,2)
imagesc(B1minusimageHR)
title 'GRE TR3000 FA90 Head Receive'


subplot(2,3,3)
imagesc(B1minusimageBR./B1minusimageHR)
caxis([0.99 1.01])
title 'Body Receive / Head Receive'

subplot(2,3,4)
imagesc((real(effAlpha).*90)./90)
caxis([0.5 1.5])
title 'effective flip angle / intended flip angle'

subplot(2,3,5)
imagesc(real(alpha))
colorbar
caxis([ 0.5 1.5])
title 'B+, Neeb2008'
subplot(2,3,6)
imagesc(real(B1minusprofileNeeb))
colorbar
caxis([0 2])
title 'B1- Profile, Neeb2008'

