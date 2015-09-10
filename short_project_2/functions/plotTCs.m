function plotTCs(FPimages, rInd, cInd, sliceNumber, offsetListNum)
% funcition to plot the times courses from fingerprinting images
TCs = zeros(numel(rInd),numel(cInd),24);
for r = 1:numel(rInd)
    for c = 1:numel(cInd)
TCs(r,c,:) = squeeze(FPimages(rInd(r),cInd(c),sliceNumber,1:24,offsetListNum));
TCs(r,c,:) = TCs(r,c,:)./TCs(r,c,1);
    end
end

figure; 
subplot 121
title 'selected time courses'
labels = cell(size(TCs,1),size(TCs,2));
for i = 1:size(TCs,1)
for j = 1:size(TCs,2)
plot(squeeze(TCs(i,j,:)),'-*')
hold on
labels{i,j} = [num2str(rInd(i)),' ',num2str(cInd(j))];
end 
end
legend (labels)
xlabel 'Image index'
ylabel 'Signal'


num2str(rInd(:));

subplot 122
title ('example image (image 1)')
imagesc(FPimages(:,:,sliceNumber,1,offsetListNum))
hold on
for r = 1:numel(rInd)
    for c = 1:numel(cInd)
plot(rInd(r),cInd(c),'r*')
text(rInd(r), cInd(c), [num2str(rInd(r)),' ',num2str(cInd(c))])
     end
end

end
