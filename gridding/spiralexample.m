load spiralexampledata
figure %JA
subplot(2,2,1) %JA
plot(kspacelocations);
subplot(2,2,2) %JA
plot(spiraldata,'.') %JA

[gdat] = gridkb(kspacelocations,spiraldata,dcf,256,1.5,2);


%figure;

% JA
subplot(2,2,3)
plot(gdat,'.')
%

%JA
subplot(2,2,4) %JA
im = fftshift(fft2(fftshift(gdat)));
im = abs(im)/500;
cmap = [0:255].'*[1 1 1] / 256;
colormap(cmap);
image(uint8(im));
colormap(cmap);
axis square;
