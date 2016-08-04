%% Use Hargeaves spiral reconstruction example to simulate ADC shift

for ADCshift = [-10:1:10]; %time points
load spiralexampledata
nInterleaves = 18;
dt = 1;

%optional simulation of an ADC shift
%sample k-space data vector 'spiraldata' at differen point
if ADCshift < 0; %start ADC before gradients
        for n = 1:size(spiraldata,2)
          spiraldata(:,n) = [zeros(abs(ADCshift)*dt,1);spiraldata(1:(end-abs(ADCshift)*dt),n)];
        end
end
if ADCshift > 0; %start ADC after gradients
        for n = 1:size(spiraldata,2)
     %       spiraldata(1+(abs(ADCshift)*dt):end,n)
                %,zeros(1,abs(ADCshift)*dt)]
            spiraldata(:,n) = [spiraldata(1+(abs(ADCshift)*dt):end,n);zeros(abs(ADCshift)*dt,1)];
        end
end


figure %JA
subplot(2,2,1) %JA
plot(kspacelocations);
title 'expected k-space locations'
subplot(2,2,2) %JA
spiraldata = spiraldata(:,1:nInterleaves);
plot(spiraldata,'.') %JA
title 'k-space data'


% grid the k-space
[gdat] = gridkb(kspacelocations,spiraldata,dcf,256,1.5,2);

%figure;

% JA
subplot(2,2,3)
plot(gdat,'.')
title 'gridded k-space data'
%

pause(0.2);
%JA
subplot(2,2,4) %JA
im = fftshift(fft2(fftshift(gdat)));
im = abs(im)/500;
cmap = [0:255].'*[1 1 1] / 256;
colormap(cmap);
image(uint8(im));
%colormap(cmap);
axis square;
title(['reconstructed image, ADC timepoint shift:',num2str(ADCshift)])

end