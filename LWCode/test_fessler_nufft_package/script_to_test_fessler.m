
ADCshifts = [0:10:110];
nChannel = 1;
for ADCshift = ADCshifts(5)
    ADCshift
filename = (['/Users/jallen/Documents/DPhil/data/TWIX/20160415/meas_MID7',num2str(28+(ADCshifts(1+(ADCshift*0.1))*0.1)),'_JA_IR_bSSFP_fp_tr10_fa10_noinv_shift',num2str(ADCshifts(1+(ADCshift*0.1))),'_FID',num2str(9497+(ADCshifts(1+(ADCshift*0.1))*0.1))]);
twix_obj = mapVBVD(filename);
twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory 
image_data = twix_obj.image{''};
kimage = squeeze(image_data(:,nChannel,:));
end
%%
load sampling_locations
im=phantom;
[ vec_of_samples]=non_uniform_fessler_fft2c(im,file);



%% All interleaves
for n = 1:48
leaf_values(((n*nSamplePts)-(nSamplePts-1)):(n*nSamplePts),1) = kimage(:,n);
spiralkspace_locations(1,(n*nSamplePts-(nSamplePts-1)):(n*nSamplePts))=spiralkspace(n,1:448);
end


vec_of_samples = 0.5*spiralkspace_locations/max(max(spiralkspace_locations));
fft_values=[real(leaf_values) imag(leaf_values) vec_of_samples'];
[ reconstructed_image ] = non_uniform_fessler_ifft2c(fft_values,size(im,1));

figure, imagesc(abs(reconstructed_image))