function [similarity, bestMatch] = calcSimilarity(data, signalDictionary, sliceNumber)
% Jack Allen.
% University of Oxford.
% jack.allen@jesus.ox.ac.uk
%
% Function to calculate the similarity of a time course with each entry in a dictionary of time courses.
%
% data(1:nRow, nColumn, sliceNumber, 1:nimageTimePts)
% similarity(1:nDataRow, 1:nDataColumn, 1:nSignalDictionaryRow, 1:nSignalDictionaryColumn)

%% simularity measure

similarity = zeros(size(data,1), size(data,2), size(signalDictionary,1), size(signalDictionary,2));
maxSimilarityScore = zeros(size(data,1), size(data,2));

for data_i = 1 : size(data,1)
    for data_j = 1 : size(data,2)
        
        for i = 1 : size(signalDictionary,1)
            for j = 1 : size(signalDictionary,2)
                similarity(data_i,data_j,i,j) = dot(squeeze(signalDictionary(i,j,:)),squeeze(data(data_i,data_j, sliceNumber, :))/(norm(squeeze(signalDictionary(i,j,:)))*norm(squeeze(data(data_i,data_j, sliceNumber,:)))));
            end
        end
        
        maxSimilarityScore(data_i,data_j) = max(max(squeeze(similarity(data_i,data_j,:,:))))
        [bestT1ind, bestT2ind] = find(squeeze(similarity(data_i,data_j,:,:)) >= maxSimilarityScore(data_i,data_j))
        matchedT1 = dictionaryParams(1, bestT1ind)
        matchedT2 = dictionaryParams(2, bestT2ind)
        
        for nbestT1ind = 1 : numel(bestT1ind)
            for nbestT2ind = 1 : numel(bestT2ind)
                bestMatch(data_i, data_j, :, nbestT1ind, nbestT2ind) = squeeze(signalDictionary(bestT1ind(nbestT1ind), bestT2ind(nbestT2ind), :))  ;
            end
        end
        
    end
    
    
    
    
    
    
    
    
end