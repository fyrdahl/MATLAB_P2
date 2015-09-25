function [similarity, matchedT1, matchedT2, matchedFAdev, M0_mean, M0_stdPC, M0,  scales, bestMatch, ind1, ind2] = calcSimilarity(data, sd, sliceNumber, dictionaryParams, workingdir)
% Jack Allen.
% University of Oxford.
% jack.allen@jesus.ox.ac.uk
%
% Function to calculate the similarity of a time course with each entry in a dictionary of time courses.
%
% data(1:nRow, nColumn, sliceNumber, 1:nimageTimePts)
% similarity(1:nDataRow, 1:nDataColumn, 1:nSignalDictionaryRow, 1:nSignalDictionaryColumn)

%% simularity measure
load([workingdir,'/MAT-files/mask.mat'])

disp('calculating similarity: started')

similarity = zeros(size(data,1), size(data,2), size(sd,1), size(sd,2), size(sd,3));
maxSimilarityScore = zeros(size(data,1), size(data,2));
matchedT1 = zeros(size(data,1), size(data,1));
matchedT2 = zeros(size(data,1), size(data,1));
matchedFAdev = zeros(size(data,1), size(data,1));
M0_mean = zeros(size(data,1), size(data,1));

tic

for data_i = 1 :size(data,1)
    for data_j = 1 : size(data,2)
        
       
        if mask(data_i,data_j) > 0
            
            for i = 1 : size(sd,1) % T1
                for j = 1 : size(sd,2)% T2
                    for k = 1 : size(sd,3) % flip angle variations due to field inhomogeneities
                        
                        
                        %                     sd = reshape(signalDictionary, [size(signalDictionary,1)*size(signalDictionary,2)*size(signalDictionary,3), size(signalDictionary,4)]);
                        %                     tc = reshape(data,[size(data,1)*size(data,2), size(data,4)]);
                        %
                        %                     similarity = dot(sd(1:end,:),tc(1:end,:))/(norm(sd)*norm(tc));
                        
                        s = squeeze(sd(i,j,k,:));
                        d = squeeze(data(data_i,data_j, sliceNumber, : ));
                        similarity(data_i,data_j,i,j,k) = dot(s,d)/(norm(s)*norm(d)) ;
                        
                    end
                    
                end
                
            end
            
            maxSimilarityScore(data_i,data_j) = max(max(max(squeeze(similarity(data_i,data_j,:,:,:)))));
            sim = squeeze(similarity(data_i,data_j,:,:,:));
            inds = find(sim == maxSimilarityScore(data_i,data_j));
            [bestT1ind, bestT2ind, bestFAdevInd] = ind2sub(size(sim),inds);
            
            matchedT1(data_i,data_j) = max(dictionaryParams(1, bestT1ind));
            matchedT2(data_i,data_j) = max(dictionaryParams(2, bestT2ind));
            matchedFAdev(data_i,data_j) = max(dictionaryParams(3, bestFAdevInd));
            
            bestMatch(data_i, data_j, :) = squeeze(sd(max(bestT1ind(:)), max(bestT2ind(:)), max(bestFAdevInd(:)) , :));
            
            for i = 1:24
            image(:,:) = data(:,:,sliceNumber,i);
            dNorm(i,1) = d(i)/max(image(:) ); %normalise data point within respect to its image
            end
            
            scales(data_i, data_j, :) = d./squeeze(bestMatch(data_i, data_j, :));
       
       % scalesNorm(data_i, data_j, :) = dNorm./squeeze(bestMatch(data_i, data_j, :));
       
            %         figure;
            %              plot(scales)
            %                 hold on
            %                 plot(d)
            %              %ylim([3000 13000])
            %              legend 'scales' 'data'
            %             pause
            
            M0_mean(data_i, data_j) =  mean(scales(data_i, data_j, :));
            M0_stdPC(data_i, data_j) = 100*(std(scales(data_i, data_j, :)))/M0_mean(data_i, data_j); %percentage
            M0(data_i, data_j) =  scales(data_i, data_j, 1);
            
           % M0_mean(data_i, data_j) =  mean(scalesNorm(data_i, data_j, :));
            %M0_stdPC(data_i, data_j) = 100*(std(scalesNorm(data_i, data_j, :)))/M0_mean(data_i, data_j); %percentage
           % M0(data_i, data_j) =  scalesNorm(data_i, data_j, 1);
 
            
            %     disp(['j percentage progress: ', num2str((data_j/size(data,2))*100)])
        end
        
        
        
        
        
    end
    
    disp(['calculating similarity: ',num2str( (data_i/size(data,1))*100) , ' percent complete'])
end



toc

end