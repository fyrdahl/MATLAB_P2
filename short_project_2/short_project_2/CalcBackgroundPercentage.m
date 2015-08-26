function[SNR signal noise] = calcSNR(images)

backgroundRoi(1,:) = [40 2 5 3] ;
signalRoi(1,:) = [24 18 5 3] ;
noise = images(backgroundRoi(1,2):(backgroundRoi(1,2)+backgroundRoi(1,4)),backgroundRoi(1,1):(backgroundRoi(1,1)+backgroundRoi(1,3)),TE(2))
noiseMean = mean(noise(:));
noiseStd = std(noise(:));

signal = images(signalRoi(1,2):(signalRoi(1,2)+signalRoi(1,4)),signalRoi(1,1):(signalRoi(1,1)+signalRoi(1,3)),TE(2));
signalMean = mean(signal(:));
signalStd = std(signal(:));

figure
imagesc(images(:,:,TE(2)))
rectangle('position',[signalRoi(1,1) signalRoi(1,2) signalRoi(1,3) signalRoi(1,4)])
rectangle('position',[backgroundRoi(1,1) backgroundRoi(1,2) backgroundRoi(1,3) backgroundRoi(1,4)])
SNR = (signal/noise)

end