
load sampling_locations
%plot sampling locations
figure, plot(file(:,1),file(:,2),'+')


im=phantom;
[ vec_of_samples]=non_uniform_fessler_fft2c(im,file);

%{
first_leaf_values=kimage(:,1);
firstspiralkspace_locations=spiralkspace(1,1:448);
vec_of_samples = 0.5*firstspiralkspace_locations/max(max(firstspiralkspace_locations));
%}



fft_values=[file vec_of_samples];
[ reconstructed_image ] = non_uniform_fessler_ifft2c(fft_values,size(im,1));

figure, 
subplot(1,2,1)
imagesc(im)
subplot(1,2,2)
imagesc(abs(reconstructed_image))

%sampling locations
figure, plot(file(:,1),file(:,2),'.')