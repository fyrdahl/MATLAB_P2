
sliceNumber = 1 % slice to be analysed
clear data

% data = FPimages(compartmentCenters(1,1),compartmentCenters(1,2),:,:);

% for r = 1:size(compartmentCenters,1)
% data(r,1,sliceNumber,:) = FPimages(compartmentCenters(r,1),compartmentCenters(r,2),sliceNumber,:);
% end

for r = 1 : size(FPimages,1)
    for c = 1 : size(FPimages,2)
data(r,c,sliceNumber,:) = FPimages(r,c,sliceNumber,:,offsetListNum);
    end
end

[similarity(:,:,offsetListNum), matchedT1(:,:,offsetListNum), matchedT2(:,:,offsetListNum), matchedFAdevInd(:,:,offsetListNum)] = calcSimilarity(data, signalDictionary(:,:,:,:,offsetListNum), sliceNumber, dictionaryParams);
%
%% visualise spread of matched T1s and T2s
% figure; hist(squeeze(matchedT1(offsetListNum,:)))
% figure;
% hist(squeeze(matchedT2(offsetListNum,:)))
% %
% figure; plot(squeeze(data(1,1,1,:)), '-*')
% hold on
% plot(squeeze(bestMatch(1,1,:))*data(1,1,1,1),'--.')
