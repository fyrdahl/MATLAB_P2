function [similarity, matchedT1, matchedT2, matchedFAdevInd, bestMatch, ind1, ind2, M0] = calcSimilarity(data, signalDictionary, sliceNumber, dictionaryParams)
% Jack Allen.
% University of Oxford.
% jack.allen@jesus.ox.ac.uk
%
% Function to calculate the similarity of a time course with each entry in a dictionary of time courses.
%
% data(1:nRow, nColumn, sliceNumber, 1:nimageTimePts)
% similarity(1:nDataRow, 1:nDataColumn, 1:nSignalDictionaryRow, 1:nSignalDictionaryColumn)

%% simularity measure
similarity = zeros(size(data,1), size(data,2), size(signalDictionary,1), size(signalDictionary,2), size(signalDictionary,3));
maxSimilarityScore = zeros(size(data,1), size(data,2));

tic
for data_i = 1 : size(data,1)
    for data_j = 1 : size(data,2)
        data_i;
        data_j;
        for i = 1 : size(signalDictionary,1) % T1
            for j = 1 : size(signalDictionary,2)% T2
                
                for k = 1 : size(signalDictionary,3) % flip angle variations due to field inhomogeneities
                    
                    s = squeeze(signalDictionary(i,j,k,:));
                    d = squeeze(data(data_i,data_j, sliceNumber, 1:numel(s) ));
                    
%                     similarity(data_i,data_j,i,j,k) = dot(s,d);
                    similarity(data_i,data_j,i,j,k) = dot(s,d)/(norm(s)*norm(d)) ;
                    
                end
                
            end
            
        end
        
                maxSimilarityScore(data_i,data_j) = max(max(max(squeeze(similarity(data_i,data_j,:,:,:)))));
                s = squeeze(similarity(data_i,data_j,:,:,:));
%                 [ind1, ind2, ind3] = find(s == maxSimilarityScore(data_i,data_j))
                inds = find(s == maxSimilarityScore(data_i,data_j));
               [bestT1ind bestT2ind bestFAdevInd] = ind2sub(size(s),inds)
          
        
%                 if mod(ind2, size(signalDictionary,3)) == 0
%                     bestT2ind = 11;
%                 else
%                     bestT2ind = mod(ind2, size(signalDictionary,3));
%                 end
%                 
%                 if ( ind3 / size(signalDictionary,3) )  < 1
%                     bestFA1devInd = 1;
%                 else
%                     bestFA1devInd = ceil(ind3./size(signalDictionary,3) );
%                 end
        
        %         matchedT1 = zeros(size(data,1), size(data,2), size(bestT1ind,1));
        %         matchedT2 = zeros(size(data,1), size(data,2), size(bestT2ind,1));
        %         matchedFA1devInd = zeros(size(data,1), size(data,2), size(bestFA1devInd,1));
        %         bestMatch = zeros(numel(bestT1ind), numel(bestT2ind), size(signalDictionary,4));
        
                matchedT1(data_i,data_j,:) = dictionaryParams(1, bestT1ind);
                matchedT2(data_i,data_j,:) = dictionaryParams(2, bestT2ind);
                matchedFAdevInd(data_i,data_j,:) = dictionaryParams(3, bestFAdevInd);
        
               
                bestMatch(data_i, data_j, :) = squeeze(signalDictionary(bestT1ind, bestT2ind, bestFAdevInd , :))  ;
        
        %         M0(data_i, data_j,:) =  squeeze(bestMatch(data_i, data_j, :, nbestT1ind, nbestT2ind, nbestFA1devInd));
  
   
    disp(['j percentage progress: ', num2str((data_j/size(data,2))*100)])
   
    end
    
    disp(['calculating similarity: ',num2str( (data_i/size(data,1))*100) , ' percent complete'])
end
toc

end