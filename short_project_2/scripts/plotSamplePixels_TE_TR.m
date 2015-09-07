figure;
imagesc(TEimages(:,:,TE(2)))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1 :plotNumCompartments
    plot(compartmentCenters(i,2,1),compartmentCenters(i,1,1),'*')
    text(compartmentCenters(i,2,1),compartmentCenters(i,1,1), compartmentLabels(i) )
end

figure;
imagesc(TIimages(:,:,TI(2)))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1:plotNumCompartments
    plot(compartmentCenters(i,2,2),compartmentCenters(i,1,2),'*')
    text(compartmentCenters(i,2,2),compartmentCenters(i,1,2), compartmentLabels(i) )
end