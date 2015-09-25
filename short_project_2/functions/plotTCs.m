function plotTCs(FPimages,data, Ind, sliceNumber, offsetListNum)


% funcition to plot the times courses from fingerprinting images
TCs = zeros(size(Ind,1),24);
for r = 1:size(Ind,1)
    
    TCs(r,:) = squeeze(FPimages(Ind(r,1),Ind(r,2),sliceNumber,1:24,offsetListNum));
    TCs(r,:) = TCs(r,:)./TCs(r,1);
    
end

figure;
%subplot 121
title 'selected time courses'
labels = cell(size(Ind,1));
for r = 1:size(Ind,1)
    
    plot(squeeze(TCs(r,:)),'-*')
    hold on
    plot(data(Ind(r,1),Ind(r,2)))
    [num2str(Ind(r,1)),' ',num2str(Ind(r,2))]
   % labels{r,:} = [num2str(Ind(r,1)),' ',num2str(Ind(r,2))];
end

%legend (labels)
xlabel 'Image index'
ylabel 'Signal'


num2str(Ind(:));

figure
%subplot 122
title ('example image (image 1)')
imagesc(FPimages(:,:,sliceNumber,1,offsetListNum))
hold on
for r = 1:size(Ind,1)
    
    plot(Ind(r,:),'r*')
    text(Ind(r,1), Ind(r,2), [num2str(Ind(r,1)),' ',num2str(Ind(r,2))])
    
end

end
