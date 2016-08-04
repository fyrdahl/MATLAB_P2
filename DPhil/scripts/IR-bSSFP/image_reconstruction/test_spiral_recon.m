run('make_spirals')
%%
for n = 1:size(spiralkspace,1)
samplelocations(1,(n-1)*nSamplePts+1:n*nSamplePts) = spiralkspace(n,1:nSamplePts);
end
samplelocations = complex(0.5*real(samplelocations)/max(real(samplelocations)),0.5*imag(samplelocations)/max(imag(samplelocations))) ;

image = phantom;
fftimage = fft(image)

imagesc(ifft(fftimage))
%generate test kspace values
rng(1);
testkspacevalues = randn(1,size(samplelocations,2));
%%
MLfile = [real(samplelocations); imag(samplelocations)]';
MLvec_of_samples = leaf_values;