function [matchedT1, matchedT2, matchedFA1devInd, bestT1ind, bestT2ind, M0] = calcSimilarity(data, signalDictionary, sliceNumber, dictionaryParams)
% Jack Allen.
% University of Oxford.
% jack.allen@jesus.ox.ac.uk
%
% Function to calculate the similarity of a time course with each entry in a dictionary of time courses.
%
% data(1:nRow, nColumn, sliceNumber, 1:nimageTimePts)
% similarity(1:nDataRow, 1:nDataColumn, 1:nSignalDictionaryRow, 1:nSignalDictionaryColumn)

%% simularity measure
signalDictionary;
similarity = zeros(size(data,1), size(data,2), size(signalDictionary,1), size(signalDictionary,2), size(signalDictionary,3));
maxSimilarityScore = zeros(size(data,1), size(data,2));

tic
for data_i = 1 : size(data,1)
    for data_j = 1 : size(data,2)
        
        for i = 1 : size(signalDictionary,1) % T1
            for j = 1 : size(signalDictionary,2)% T2
                for k = 1 : size(signalDictionary,3) % flip angle variations due to field inhomogeneities
                    A = squeeze(signalDictionary(i,j,k,:));  
                    B = squeeze(data(data_i,data_j, sliceNumber, 1:numel(A) ));
                    similarity(data_i,data_j,i,j,k) = dot(A,B)/(norm(squeeze(signalDictionary(i,j,k,:)))*norm(squeeze(data(data_i,data_j, sliceNumber,:))));
                end
            end
        end
        
        maxSimilarityScore(data_i,data_j) = max(max(max(squeeze(similarity(data_i,data_j,:,:,:)))))
%         s = squeeze(similarity(data_i,data_j,:,:,:))
%         [bestT1ind, b] = find(s == maxSimilarityScore(data_i,data_j));
%         
%         if mod(b, size(signalDictionary,3)) == 0            
%             bestT2ind = 11;
%         else
%             bestT2ind = mod(b, size(signalDictionary,3));
%         end
%         
%        
%         if ( b / size(signalDictionary,3) )  < 1
%         bestFA1devInd = 1
%         else
%         bestFA1devInd = ceil(b./size(signalDictionary,3) )
%         end
%             
%         matchedT1(data_i,data_j,:) = dictionaryParams(1, bestT1ind);
%         matchedT2(data_i,data_j,:) = dictionaryParams(2, bestT2ind);
%         matchedFA1devInd(data_i,data_j,:) = dictionaryParams(3, bestFA1devInd);
%         
%         for nbestT1ind = 1 : numel(bestT1ind)
%             for nbestT2ind = 1 : numel(bestT2ind)
%                 for nbestFA1devInd = 1 : numel(bestFA1devInd)
%                     bestMatch(data_i, data_j, :, nbestT1ind, nbestT2ind, nbestFA1devInd) = squeeze(signalDictionary(bestT1ind(nbestT1ind), bestT2ind(nbestT2ind), bestFA1devInd(nbestFA1devInd) , :))  ;
%                 end
%             end
%         end
%         
%         
%         M0(data_i, data_j,:) =  squeeze(bestMatch(data_i, data_j, :, nbestT1ind, nbestT2ind, nbestFA1devInd))
    end
    
    disp(['calculating similarity: ',num2str( (data_i/size(data,1))*100) , ' percent complete'])
end
toc

end