function[SNR,bgStd] = calcSNR(images,t,bgROI,sigROI,figureFlag)

backgroundRoi(1,:) = bgROI;
signalRoi(1,:) = sigROI;
background = images(backgroundRoi(1,2):(backgroundRoi(1,2)+backgroundRoi(1,4)),backgroundRoi(1,1):(backgroundRoi(1,1)+backgroundRoi(1,3)),t(2));
bgStd = std(background(:));

signal = images(signalRoi(1,2):(signalRoi(1,2)+signalRoi(1,4)),signalRoi(1,1):(signalRoi(1,1)+signalRoi(1,3)),t(2));
signalMean = mean(signal(:));
SNR = 0.65*(signalMean/bgStd);

if figureFlag == 1
figure
imagesc(images(:,:,t(2)))
rectangle('position',[signalRoi(1,1) signalRoi(1,2) signalRoi(1,3) signalRoi(1,4)])
rectangle('position',[backgroundRoi(1,1) backgroundRoi(1,2) backgroundRoi(1,3) backgroundRoi(1,4)])
end
end