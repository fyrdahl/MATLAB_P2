function [similarity, matchedT1, matchedT2, matchedFAdev, M0_mean, M0_stdPC, M0, M0fit_grad,  scales, bestMatch, el, ind1, ind2] = calcSimilarity(data, signalDictionary, sliceNumber, dictionaryParams, workingdir)
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

data = reshape(data,[64*64, 24]);
mask = reshape(mask,[64*64, 1]);
sd = reshape(signalDictionary, [ size(signalDictionary,1)*size(signalDictionary,2)*size(signalDictionary,3), 24]);

similarity = zeros(size(data,1), size(sd,1) );
maxSimilarityScore = zeros(size(data,1),1);
matchedT1 = zeros(size(data,1),1);
matchedT2 = zeros(size(data,1),1);
matchedFAdev = zeros(size(data,1),1);
M0_mean = zeros(size(data,1),1);
scales = zeros(size(data,1),size(data,2));
bestMatch = zeros(size(data,1),size(data,2));
M0fit_grad = zeros(size(data,1),1);

M0model = @(a,x) a*x;

tic

pp = parpool(4)
for data_i = 1 :size(data,1)
    %for data_j = 1 : size(data,2)
    
    
    if mask(data_i,1) > 0
        
        parfor sd_i = 1 : size(sd,1) % T1
            
            
            %                     sd = reshape(signalDictionary, [size(signalDictionary,1)*size(signalDictionary,2)*size(signalDictionary,3), size(signalDictionary,4)]);
            %                     tc = reshape(data,[size(data,1)*size(data,2), size(data,4)]);
            %
            %                     similarity = dot(sd(1:end,:),tc(1:end,:))/(norm(sd)*norm(tc));
            
            
            similarity(data_i,sd_i) = dot(sd(sd_i,:),data(data_i,:))/(norm(sd(sd_i,:))*norm(data(data_i,:))) ;
            
            
            
        end
        
        maxSimilarityScore(data_i) = max(similarity(data_i,:));
        sim = similarity(data_i,:);
        sim = reshape(sim,[size(signalDictionary,1),size(signalDictionary,2),size(signalDictionary,3)]);
        inds = find(sim == maxSimilarityScore(data_i));
       
        [bestT1ind, bestT2ind, bestFAdevInd] = ind2sub(size(sim),inds);
        
      matchedT1(data_i) = max(dictionaryParams(1, bestT1ind));
      matchedT2(data_i) = max(dictionaryParams(2, bestT2ind));
      matchedFAdev(data_i) = max(dictionaryParams(3, bestFAdevInd));
        
      
       bestMatch(data_i, :) = squeeze(signalDictionary(max(bestT1ind(:)), max(bestT2ind(:)), max(bestFAdevInd(:)) , :));
       
       scales(data_i, :) = data(data_i,:)./squeeze(bestMatch(data_i, :));
        
       M0fit = fit(squeeze(bestMatch(data_i, :))', data(data_i,:)',M0model,'Upper',[6000],'Lower',[0],'StartPoint',[1000] );
       
       M0fit_grad(data_i,1) = M0fit.a;
  
        % scalesNorm(data_i, data_j, :) = dNorm./squeeze(bestMatch(data_i, data_j, :));
        
        %         figure;
        %              plot(scales)
        %                 hold on
        %                 plot(d)
        %              %ylim([3000 13000])
        %              legend 'scales' 'data'
        %             pause
     
        M0_mean(data_i) =  mean(scales(data_i, :));
        M0_stdPC(data_i) = 100*(std(scales(data_i, :)))/M0_mean(data_i); %percentage
        M0(data_i) =  scales(data_i, 1);
        
        % M0_mean(data_i, data_j) =  mean(scalesNorm(data_i, data_j, :));
        %M0_stdPC(data_i, data_j) = 100*(std(scalesNorm(data_i, data_j, :)))/M0_mean(data_i, data_j); %percentage
        % M0(data_i, data_j) =  scalesNorm(data_i, data_j, 1);
        
        
        %     disp(['j percentage progress: ', num2str((data_j/size(data,2))*100)])
    end
    
    
    
    
    
    
    
    disp(['calculating similarity: ',num2str( (data_i/size(data,1))*100) , ' percent complete'])
    
end
matchedT1 = reshape(matchedT1, [sqrt(size(matchedT1)), sqrt(size(matchedT1))]);
matchedT2 = reshape(matchedT2, [sqrt(size(matchedT2)), sqrt(size(matchedT2))]);
matchedFAdev = reshape(matchedFAdev, [sqrt(size(matchedFAdev)), sqrt(size(matchedFAdev))]);
M0_mean = reshape(M0_mean, [sqrt(size(M0_mean)), sqrt(size(M0_mean))]);
M0fit_grad = reshape(M0fit_grad, [sqrt(size(data,1)),sqrt(size(data,1))]);
delete(pp) %shutdown parpool
el = toc

end