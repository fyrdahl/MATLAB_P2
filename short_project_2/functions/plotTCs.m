function plotTCs(data)
%% Jack Allen
% jack.allen@jesus.ox.ac.uk
%
% function to plot the times courses from specific coordinates from images
% acquired at different time points.
%
%% Select specific points
fig = figure;
subplot 121
title ('example image (image 1)')
testTC = data(size(data,1)/2,size(data,2)/2,:);
testInd = find(testTC);
imagesc(data(:,:,testInd(1)))
[x, y] = getpts(fig);
x = round(x);
y = round(y);
nTimePts = size(data,3);
TCs = zeros(size(x,1),nTimePts);
for i = 1:size(x,1)   
    TCs(i,:) = squeeze(data(x(i),y(i),:));
   % TCs(i,:) = TCs(i,:)./TCs(i,1);
end
%% overlay labels of selected points
title ('example image (image 1)')
imagesc(data(:,:,testInd(1)))
hold on
for i = 1:size(x,1)  
    plot(x(i),y(i),'r*')
    text(x(i), y(i), [num2str(x(i)),' ',num2str(y(i))])
   labels(i,:) = [x(i), y(i)];
end
%% Plot the timecourses for the selected points
subplot 122
title 'selected time courses'
for i = 1:size(x,1)  
    plot(TCs(i,:),'-*')
    hold on
end
legend([num2str(labels)])
xlabel 'Image Index'
ylabel 'Normalised Signal'

end
