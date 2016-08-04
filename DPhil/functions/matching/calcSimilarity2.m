function [maxSimilarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime] = calcSimilarity2(data, sd, dictionaryParams)
% calcSimilarity
% Description: A function to calculate the similarity of a time course with
% each entry in a dictionary of time courses and find the best match.
%
% Author: Jack Allen.
% Institution: University of Oxford.
% Contact: jack.allen@jesus.ox.ac.uk
%
%% simularity measure
disp('calculating similarity [calSimilarity2]: started')

%mask to exclude time course outside of ROI
%loadMask.Mask = load([savingdir,'/MAT-files/mask.mat']);
%mask = reshape(loadMask.mask,[64*64, 1]);

nPts = size(data,3);
data = reshape(data,size(data,1)*size(data,2), nPts);

%initialise maps
tempTCsimilarity = zeros(size(sd,1),1);
maxSimilarity = zeros(size(data,1),1);
bestMatch = zeros(size(data,1),nPts);
matchedT1 = zeros(size(data,1),1);
matchedT2 = zeros(size(data,1),1);
matchedFAdev = zeros(size(data,1),1);
M0fit_grad = zeros(size(data,1),1);

M0model = @(a,x) a*x;

tic % start timing how long matching takes

disp('calculating similarity...')
for data_i = 1:size(data,1) %for all pixels
    
    %compare the TC with each dictionary entry
    %if mask(data_i,1) > 0
    for sd_i = 1:size(sd,1) %for all entries in the dictionary
        tempTCsimilarity(sd_i) = dot(sd(sd_i,:), data(data_i,:))/(norm(sd(sd_i,:))*norm(data(data_i,:)));
    end
    
    %find highest similarity score for the TC of interest
    [maxSimilarity(data_i)]  = max(tempTCsimilarity(:));
    
    %Find the parameters associated with the best match
    tempTCsimilarity =  reshape(tempTCsimilarity,[size(signalDictionary,1),size(signalDictionary,2),size(signalDictionary,3)]);
    
    [bestT1ind, bestT2ind, bestFAdevInd] = ind2sub(size(tempTCsimilarity), find(tempTCsimilarity == maxSimilarity(data_i)));
    
    %Assign the parameters  to the timecourse of interest
    matchedT1(data_i,1) = dictionaryParams(1, bestT1ind);
    matchedT2(data_i,1) = dictionaryParams(2, bestT2ind);
    matchedFAdev(data_i,1) = dictionaryParams(3, bestFAdevInd);
    bestMatch(data_i, :) = squeeze(signalDictionary(bestT1ind, bestT2ind, bestFAdevInd, 1:nTimeCoursePts));
    M0fit = fit(squeeze(bestMatch(data_i, :))', data(data_i,:)',M0model,'Upper',6000,'Lower',0,'StartPoint',500 );
    M0fit_grad(data_i,1) = M0fit.a;
    
    disp(['calculating similarity: ',num2str( (data_i/size(data,1))*100) , ' percent complete'])
    
end
%reshape the assigned parameter maps into 2D form
maxSimilarity = reshape(maxSimilarity,[sqrt(size(maxSimilarity,1)),sqrt(size(maxSimilarity,1))]);
matchedT1 = reshape(matchedT1, [sqrt(size(matchedT1,1)), sqrt(size(matchedT1,1))]);
matchedT2 = reshape(matchedT2, [sqrt(size(matchedT2,1)), sqrt(size(matchedT2,1))]);
matchedFAdev = reshape(matchedFAdev, [sqrt(size(matchedFAdev,1)), sqrt(size(matchedFAdev,1))]);
M0fit_grad = reshape(M0fit_grad, [sqrt(size(data,1)),sqrt(size(data,1))]);
bestMatch = reshape(bestMatch, [sqrt(size(bestMatch,1)), sqrt(size(bestMatch,1)),nTimeCoursePts]);


matchTime = toc;

end