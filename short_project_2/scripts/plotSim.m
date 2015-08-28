
figure; hold on
% title (['fingerprint offset list ', num2str(offsetListNum)])
ySim = simImageMtransverse(:)./simImageMtransverse(1);
ySim = Mxy(:)./Mxy(1);
% plot(ySim,'x','MarkerSize',20)
plot(ySim,'x','MarkerSize',20)
y = zeros(3,(size(FPimages,4)/2));

sliceNumber = 1;

%mean of rectangle covering most of the phantom
for i = 1 : size((FPimages(25:35,25:35,sliceNumber,:)) ,4)/2
    ROI = (FPimages(25:35,25:35,sliceNumber,i));
    ROIimage(i,:) = ROI(:);
    ROIimageMean(i) = mean(ROIimage(i,:));
    ROIimageStd(i) = std(ROIimage(i,:));
end
% % rectangle('position',[25 25 10 10])
% %
normROIimageMean = squeeze(ROIimageMean/ROIimageMean(1));
yImageROI = normROIimageMean;
% plot(yImageROI(1:24),'^')
errorbar(yImageROI(1:24),ROIimageStd/ROIimageMean(1),'^');

%signal from each compartment
title (['Phantom: ',phantomName,', Offset list: ',num2str(offsetListNum),', compartment center coords list: ',num2str(compartmentCentersList)]);
for n = 1:plotNumCompartments
    for i = 1:(size(FPimages,4)/2)
        y(n,i) = squeeze(FPimages(compartmentCenters(n,1,compartmentCentersList),compartmentCenters(n,2,compartmentCentersList),sliceNumber,i));
    end
    normStdBG = (std(background(:)))/y(n,1);
    y(n,:) = y(n,:)/y(n,1);
    % residuals = y(n,:) - ySim;
%     plot(y(n,1:(size(FPimages,4)/2)),'*')
    
    errorbar(y(n,:),repmat(normStdBG,1,24),'*' );
end
legend 'Simulated Signal' 'ROI mean and std' 'pixel with BG std' 'pixel with BG std' 
xlabel 'TE indices'
ylabel 'signal (normalised to first measurement)'


%  figure; plot(residuals,'+')