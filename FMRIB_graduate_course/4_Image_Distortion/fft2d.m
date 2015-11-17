function out = fft2d(in)

out = fftshift(fft(fft(ifftshift( in ),[],1),[],2));