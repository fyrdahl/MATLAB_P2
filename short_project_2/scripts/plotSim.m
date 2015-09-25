
simCom = figure; hold on
% title (['fingerprint offset list ', num2str(offsetListNum)])
ySim = Mxy(:)./Mxy(1);
plot(ySim,'x','MarkerSize',20)
y = zeros(3,(size(FPimages,4)/2));

sliceNumber = 1;

%% mean of rectangle covering most of the phantom
% for i = 1 : size((FPimages(25:35,25:35,sliceNumber,:)) ,4)/2
%     ROI = (FPimages(25:35,25:35,sliceNumber,i));
%     ROIimage(i,:) = ROI(:);
%     ROIimageMean(i) = mean(ROIimage(i,:));
%     ROIimageStd(i) = std(ROIimage(i,:));
% end
% % rectangle('position',[25 25 10 10])
% %
% normROIimageMean = squeeze(ROIimageMean/ROIimageMean(1));
% yImageROI = normROIimageMean;
% % plot(yImageROI(1:24),'^')
% errorbar(yImageROI(1:24),ROIimageStd/ROIimageMean(1),'^');

%% signal from each compartment
%title (['Phantom: ',phantomName,', Offset list: ',num2str(offsetListNum),', compartment center coords list: ',num2str(compartmentCentersList)]);
for n = 1:plotNumCompartments
    for i = 1:(size(FPimages,4)/2)
        y(n,i) = squeeze(FPimages(compartmentCenters(n,1,compartmentCentersList),compartmentCenters(n,2,compartmentCentersList),sliceNumber,i,offsetListNum));
    end
    normStdBG = (std(background(:)))/y(n,1);
    y(n,:) = y(n,:)/y(n,1);
    % residuals = y(n,:) - ySim;
    plot(y(n,1:(size(FPimages,4)/2)),'.')
    
   % errorbar(y(n,:),repmat(normStdBG,1,24),'.' );
   
end
ylim([0 1.1])
legend ({'Simulated Signal', 'Sample Pixel 1','Sample Pixel 2','Sample Pixel 3','Sample Pixel 4','Sample Pixel 5','Sample Pixel 6'},'Position',[0.6, 0.7, 0.1,0.1])
xlabel 'TE Index'
ylabel 'Normalised Signal'
savefig([workingdir,'/figures/compareSimwithData_Phantom_',phantomName,'__Offset_list_',num2str(offsetListNum),'_compartmentcentercoordslist:',num2str(compartmentCentersList),'.fig'])

matlab2tikz('figurehandle',simCom,'filename',[filename,phantomName,'simCom',num2str(offsetListNum)],'height', '\figureheight', 'width', '\figurewidth')
    

%% figure; plot(residuals,'+')