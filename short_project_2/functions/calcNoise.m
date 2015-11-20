function[bgMean,bgStd,SNR] = calcNoise(images,t,bgROI,sigROI,figureFlag)
% calcNoise
% 
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
%
%

backgroundRoi(1,:) = bgROI;
signalRoi(1,:) = sigROI;
background = images(backgroundRoi(1,2):(backgroundRoi(1,2)+backgroundRoi(1,4)),backgroundRoi(1,1):(backgroundRoi(1,1)+backgroundRoi(1,3)),t);
bgStd = std(background(:));
bgMean = mean(background(:));

signal = images(signalRoi(1,2):(signalRoi(1,2)+signalRoi(1,4)),signalRoi(1,1):(signalRoi(1,1)+signalRoi(1,3)),t);
signalMean = mean(signal(:));
SNR = 0.65*(signalMean/bgStd);

if figureFlag == 1
figure
imagesc(images(:,:,t))
rectangle('position',[signalRoi(1,1) signalRoi(1,2) signalRoi(1,3) signalRoi(1,4)])
rectangle('position',[backgroundRoi(1,1) backgroundRoi(1,2) backgroundRoi(1,3) backgroundRoi(1,4)])
end
end