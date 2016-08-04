load sampling_locations
im=phantom;
[ vec_of_samples]=non_uniform_fessler_fft2c(im,file);
fft_values=[file vec_of_samples];
 [ reconstructed_image ] = non_uniform_fessler_ifft2c(fft_values,size(im,1));