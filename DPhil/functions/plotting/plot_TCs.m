function plot_TCs(data,varargin)
% function plot_TCs(data,varargin)
%
% DATA: 3D matrix
% PHANTOMNAME: string
% OFFSETLISTNUM: number
% VARARGIN{1}: specify whether or not the data is fingerprinting data
% (set as 'FP' if so)
% VARARGIN{2}: Best matches
%
% Function to plot the time courses from selected pixels, with the option
% of also plotting the best match from a dictionary if the images were
% obtained via a fingerprinting method.
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>

% Select specific points
fig = figure;
subplot 121
title ('example image (image 1)')
testTC = data(size(data,1)/2,size(data,2)/2,:);
testInd = find(testTC);
imagesc(data(:,:,testInd(1)))
[row, col] = getpts(fig);
row = round(row);
col = round(col);
nTimePts = size(data,3);
dataTCs = zeros(size(row,1),nTimePts);

set(gca,'FontSize',24)
for i = 1:size(row,1)
    dataTCs(i,:) = squeeze(data(row(i),col(i),:));
    x(i,:) = find(dataTCs(i,:)>0);
if ~isempty(varargin)
    dataTCs(i,:) = dataTCs(i,:)./dataTCs(i,x(i,1));
end
end

set(gca,'FontSize',24)
%% overlay labels of selected points
title ('example image (image 1)')
imagesc(data(:,:,testInd(1)))
hold on
for i = 1:size(row,1)
    plot(row(i),col(i),'r*')
     set(gca,'FontSize',24)
    text(row(i), col(i), [num2str(row(i)),' ',num2str(col(i))])
    labels{1,i} = [num2str(row(i)),' ',num2str(col(i))];

end
%% Plot the timecourses for the selected points
subplot 122
set(gca,'FontSize',24)
title 'selected time courses'
for i = 1:size(row,1)
    plot(x(i,:),dataTCs(i,x(i,:)),'-o')
    set(gca,'FontSize',24)
    hold on
end
legend(labels)
xlabel('Image Index')
ylabel('Normalised Intensity')



if ~isempty(varargin)
% varargin{3} = best matches
bestMatch = varargin{2};
bestMatchTCs = zeros(size(row,1),nTimePts);
for i = 1:size(row,1)
    set(gca,'FontSize',24)
    bestMatchTCs(i,:) = squeeze(bestMatch(row(i),col(i),:));
    bestMatchTCs(i,:) = bestMatchTCs(i,:)./bestMatchTCs(i,1);
    bestMatchLabels{1,i} = [labels{1,i},' best match'];
    plot(bestMatchTCs(1,:),'--*');
    set(gca,'FontSize',24,'LineWidth',2)
end
    legend([labels,bestMatchLabels])
    set(gca,'FontSize',24)
end

xlabel 'Image Index'
ylabel 'Normalised Signal'

end
