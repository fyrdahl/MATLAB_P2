
figure;
imagesc(FPimages(:,:,1))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1:plotNumCompartments
    plot(compartmentCenters(i,1,compartmentCentersList),compartmentCenters(i,2,compartmentCentersList),'*')
    text(compartmentCenters(i,1,compartmentCentersList),compartmentCenters(i,2,compartmentCentersList), compartmentLabels(i) )
end