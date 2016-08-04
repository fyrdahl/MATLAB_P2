%% 1. Plot channel data
subplot 221
imagesc(phaseImg(:,:,1))
title 'Phase Image - Ch 1'
subplot 222
imagesc(phaseImg(:,:,2))
title 'Phase Image - Ch 2'
subplot 223
imagesc(phaseImg(:,:,3))
title 'Phase Image - Ch 3'
subplot 224
imagesc(phaseImg(:,:,4))
title 'Phase Image - Ch 4'

figure
subplot 221
imagesc(magImg(:,:,1))
title 'Magnitude - Ch 1'
subplot 222
imagesc(magImg(:,:,2))
title 'Magnitude - Ch 2'
subplot 223
imagesc(magImg(:,:,3))
title 'Magnitude - Ch 3'
subplot 224
imagesc(magImg(:,:,4))
title 'Magnitude - Ch 4'

%% 2.1 Reconstruct channel data
% read in magnitude and phase images for all channels
phaseFile = '/Users/jallen/Documents/DPhil/raw_data/NifTI/20160304_MR_PHYSICS_435/1_008_JA_IR_bSSFP_fp_matrix_20160304/20160304_008_JA_IR_bSSFP_fp_matrix.nii';
[phaseImg] = read_avw(phaseFile);
phaseImg = squeeze(phaseImg(:,:,1,:));
magFile = '/Users/jallen/Documents/DPhil/raw_data/NifTI/20160304_MR_PHYSICS_435/1_006_JA_IR_bSSFP_fp_matrix_20160304/20160304_006_JA_IR_bSSFP_fp_matrix.nii';
[magImg] = read_avw(magFile);
magImg = squeeze(magImg(:,:,1,:));

combinedImage = zeros(256,256);
images = zeros(256,256,4);

for chNum = 1:4;

% check that the scanner reconstruction factor was low enough
mag = magImg(:,:,chNum);
phase = phaseImg(:,:,chNum);


maxMag = max(mag(:));
maxPhase = max(phase(:));

%phaseImg(:,:,chNum) = (360/maxPhase)*phase;
phase = ((pi)/maxPhase)*phase;

% Use the magnitude and phae information to produce a complex image
fftimage = mag.*exp(i*phase);



% do an inverse fourier transform on the complex data
image = ifft(fftimage,[],2);

% take magnitude 
images(:,:,chNum) = abs(image);

imageMin(chNum) = min(min(abs(image)));
imageMax(chNum) = max(max(abs(image)));

combinedImage(:,:) = combinedImage+abs(image);
end
%% 2.2 plot reconstruction results 
figure
subplot 221
imagesc(images(:,:,1))
title('ch1')

subplot 222
imagesc(images(:,:,2))
title('ch2')

subplot 223
imagesc(images(:,:,3))
title('ch3')

subplot 224
imagesc(images(:,:,4))
title('ch4')

figure
imagesc(combinedImage)

%% 3. Test Fourier transforms
image=imread('cameraman.tif');
figure, imshow(image), colormap gray
title('original image')

fftimage = fft2(double(image));
mag = abs(fftimage);
phase = angle(fftimage);

fftimage_custom = mag.*exp(i*phase);

%figure, imshow(abs(fftshift(fftimage_custom)),[0 100000]), colormap gray
%title('fftimage Magnitude')

ifftImg = ifft2(fftimage_custom);
imageMin = min(min(abs(ifftImg)));
imageMax = max(max(abs(ifftImg)));
figure, imshow(abs(ifftImg), [imageMin imageMax]), colormap gray
title('ifftimage  Magnitude')


