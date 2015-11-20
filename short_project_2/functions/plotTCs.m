function plotTCs(data,savingdir,phantomName,offsetListNum,FPflag)
% plotTCs: Plots of signal time courses
%
% Function to plot the time courses from selected pixels, with the option
% of also plotting the best match from the dictionary, if the images were
% obtained via the fingerprinting method.
%
% plotTCs(data,savingdir,phantomName,offsetListNum,FPflag)
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>

% Select specific points
fig = figure;
subplot 121
title ('example image (image 1)')
testTC = data(size(data,1)/2,size(data,2)/2,:);
bestMatchLoad = load([savingdir,'/MAT-files/matches/BestMatch/',phantomName,'offsetlist',num2str(offsetListNum),'paramList3bestMatch.mat'],'bestMatch');
bestMatch = bestMatchLoad.bestMatch;
testInd = find(testTC);
imagesc(data(:,:,testInd(1)))
[x, y] = getpts(fig);
x = round(x);
y = round(y);
nTimePts = size(data,3);
dataTCs = zeros(size(x,1),nTimePts);
bestMatchTCs = zeros(size(x,1),nTimePts);
for i = 1:size(x,1)   
    dataTCs(i,:) = squeeze(data(x(i),y(i),:));
    bestMatchTCs(i,:) = squeeze(bestMatch(x(i),y(i),:));
   dataTCs(i,:) = dataTCs(i,:)./dataTCs(i,1);
   bestMatchTCs(i,:) = bestMatchTCs(i,:)./bestMatchTCs(i,1);
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
    plot(dataTCs(i,:),'-')
    hold on
end
legend([num2str(labels)])

if FPflag == '1'
if size(x,1) == 1 %if only one point is selected, the best match will also be plotted.
plot(bestMatchTCs(1,:),'r-');

legend(num2str(labels),'Best Match')
end
end
xlabel 'Image Index'
ylabel 'Normalised Signal'

end
