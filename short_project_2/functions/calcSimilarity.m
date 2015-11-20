function [similarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, el] = calcSimilarity(data, nTimeCoursePts, dictionaryParams, paramRangeList,savingdir, phantomName,offsetListNum)
% Author: Jack Allen.
% Institution: University of Oxford.
% Contact: jack.allen@jesus.ox.ac.uk
%
% Description: A function to calculate the similarity of a time course with
% each entry in a dictionary of time courses and find the best match.

%% simularity measure
disp('calculating similarity: started')

load([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramRangeList),'dictionary.mat']);
signalDictionary = sd(:,:,:,:,offsetListNum);
signalDictionary = squeeze(signalDictionary);
%mask to exclude time course outside of ROI
%loadMask.Mask = load([savingdir,'/MAT-files/mask.mat']);
%mask = reshape(loadMask.mask,[64*64, 1]);

%time course of interest
data = reshape(data,[size(data,1)*size(data,2), nTimeCoursePts]);
%simulated signals
sd = reshape(signalDictionary, [ size(signalDictionary,1)*size(signalDictionary,2)*size(signalDictionary,3), nTimeCoursePts]);

%initialise maps
similarity = zeros(size(data,1), size(sd,1) );
maxSimilarityScore = zeros(size(data,1),1);
bestMatch = zeros(size(data,1),nTimeCoursePts);
matchedT1 = zeros(size(data,1),1);
matchedT2 = zeros(size(data,1),1);
matchedFAdev = zeros(size(data,1),1);
M0fit_grad = zeros(size(data,1),1);

M0model = @(a,x) a*x;

tic

pp = parpool(4);
for data_i = 1:size(data,1) %for all pixels
    disp('calculating similarity...')
    
    %if mask(data_i,1) > 0
    parfor sd_i = 1:size(sd,1) % for all entries in the dictionary
        tempData = data;
        TCsimilarity(data_i,sd_i) = dot(sd(sd_i,:), tempData(data_i,:))/(norm(sd(sd_i,:))*norm(tempData(data_i,:)));
    end
    TCsimilarity = reshape(TCsimilarity,[size(signalDictionary,1),size(signalDictionary,2),size(signalDictionary,3)]);
    
    %   find highest similarity score for TC of interest
    [maxSimilarityScore(data_i), bestSimulatedSignalIndex]  = max(TCsimilarity(:));
    %Find the parameters associated with the best match
    
    [bestT1ind, bestT2ind, bestFAdevInd] = ind2sub(size(TCsimilarity),bestSimulatedSignalIndex);
    
    %Assign the parameters to the timecourse of interest
    matchedT1(data_i,1) = dictionaryParams(1, bestT1ind);
    matchedT2(data_i,1) = dictionaryParams(2, bestT2ind);
    matchedFAdev(data_i,1) = dictionaryParams(3, bestFAdevInd);
    bestMatch(data_i, :) = squeeze(signalDictionary(bestT1ind, bestT2ind, bestFAdevInd, :));
    M0fit = fit(squeeze(bestMatch(data_i, :))', data(data_i,:)',M0model,'Upper',6000,'Lower',0,'StartPoint',500 );
    M0fit_grad(data_i,1) = M0fit.a;
    
    % end
    disp(['calculating similarity: ',num2str( (data_i/size(data,1))*100) , ' percent complete'])
    
end

%reshape the assigned parameter maps into 2D form
matchedT1 = reshape(matchedT1, [sqrt(size(matchedT1,1)), sqrt(size(matchedT1,1))]);
matchedT2 = reshape(matchedT2, [sqrt(size(matchedT2,1)), sqrt(size(matchedT2,1))]);
matchedFAdev = reshape(matchedFAdev, [sqrt(size(matchedFAdev,1)), sqrt(size(matchedFAdev,1))]);
M0fit_grad = reshape(M0fit_grad, [sqrt(size(data,1)),sqrt(size(data,1))]);
bestMatch = reshape(bestMatch, [sqrt(size(bestMatch,1)), sqrt(size(bestMatch,1)),nTimeCoursePts]);

%shutdown parpool
delete(pp)

el = toc;

save([savingdir,'/MAT-files/matches/similarity/',phantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'similarity.mat'],'similarity')
save([savingdir,'/MAT-files/matches/T1/',phantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'matchedT1.mat'],'matchedT1')
save([savingdir,'/MAT-files/matches/T2/',phantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'matchedT2.mat'],'matchedT2')
save([savingdir,'/MAT-files/matches/B1/',phantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'matchedFAdev.mat'],'matchedFAdev')
save([savingdir,'/MAT-files/matches/BestMatch/',phantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'bestMatch.mat'],'bestMatch')
save([savingdir,'/MAT-files/matches/MatchingTimes/',phantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'compileDictionaryElapsedTime.mat'],'match_time')
save([savingdir,'/MAT-files/matches/M0/',phantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'M0fit_grad.mat'],'M0fit_grad')
end