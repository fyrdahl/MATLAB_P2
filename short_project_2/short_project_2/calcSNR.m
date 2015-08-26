function[SNR signal background] = calcSNR(images,t,figureFlag)

backgroundRoi(1,:) = [40 2 7 7] ;
signalRoi(1,:) = [24 18 7 7] ;
background = images(backgroundRoi(1,2):(backgroundRoi(1,2)+backgroundRoi(1,4)),backgroundRoi(1,1):(backgroundRoi(1,1)+backgroundRoi(1,3)),t(2));
backgroundMean = mean(background(:));
backgroundStd = std(background(:));

signal = images(signalRoi(1,2):(signalRoi(1,2)+signalRoi(1,4)),signalRoi(1,1):(signalRoi(1,1)+signalRoi(1,3)),t(2));
signalMean = mean(signal(:));
signalStd = std(signal(:));

if figureFlag == 'showFigure'
figure
imagesc(images(:,:,t(2)))
rectangle('position',[signalRoi(1,1) signalRoi(1,2) signalRoi(1,3) signalRoi(1,4)])
rectangle('position',[backgroundRoi(1,1) backgroundRoi(1,2) backgroundRoi(1,3) backgroundRoi(1,4)])
SNR = (signal./background);
end
end