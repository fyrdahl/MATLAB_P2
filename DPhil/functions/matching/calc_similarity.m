<<<<<<< HEAD
function [matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime, maxSimilarity, degeneracyMap] = calc_similarity(data, dictionary, dictionaryParams)
% [maxSimilarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime] = calc_similarity(data, nTimeCoursePts, dictionaryParams, paramRangeList,savingdir, inputPhantomName,outputPhantomName,offsetListNum)
=======
function [maxSimilarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime] = calcSimilarity(data, nTimeCoursePts, dictionaryParams, paramRangeList,savingdir, inputPhantomName,outputPhantomName,offsetListNum)
% calcSimilarity
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
% Description: A function to calculate the similarity of a time course with
% each entry in a dictionary of time courses and find the best match.
%
% Author: Jack Allen.
% Institution: University of Oxford.
% Contact: jack.allen@jesus.ox.ac.uk
%
<<<<<<< HEAD
%
%% Similarity measurements for magnetic resonance fingerprinting (MRF) matching
disp('calculating similarity: started')

=======
%% simularity measure
disp('calculating similarity: started')

load([savingdir,'/MAT-files/dictionaries/',inputPhantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramRangeList),'dictionary.mat']);
size(sd)
signalDictionary = sd(:,:,:,:,offsetListNum);
signalDictionary = squeeze(signalDictionary);
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
%mask to exclude time course outside of ROI
%loadMask.Mask = load([savingdir,'/MAT-files/mask.mat']);
%mask = reshape(loadMask.mask,[64*64, 1]);

<<<<<<< HEAD
nTimeCoursePts = size(data,3);

dataDim1 = size(data,1);
dataDim2 = size(data,2);
%time course of interest
data = reshape(data,[size(data,1)*size(data,2), nTimeCoursePts]);

%% Initialise maps
maxSimilarity = zeros(size(data,1),1);
degeneracyMap = zeros(size(data,1),1);
=======
%time course of interest
data = reshape(data,[size(data,1)*size(data,2), nTimeCoursePts]);
%simulated signals
sd = reshape(signalDictionary(:,:,:,1:nTimeCoursePts), [ size(signalDictionary,1)*size(signalDictionary,2)*size(signalDictionary,3), nTimeCoursePts]);

%initialise maps
tempTCsimilarity = zeros(size(sd,1),1);
maxSimilarity = zeros(size(data,1),1);
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
bestMatch = zeros(size(data,1),nTimeCoursePts);
matchedT1 = zeros(size(data,1),1);
matchedT2 = zeros(size(data,1),1);
matchedFAdev = zeros(size(data,1),1);
M0fit_grad = zeros(size(data,1),1);

<<<<<<< HEAD
M0model = @(a,x)a*x; % for calculating M0

nParam1 = numel(dictionaryParams(1,dictionaryParams(1,:)>0)); % number of different T1
nParam2 = numel(dictionaryParams(1,dictionaryParams(2,:)>0)); % number of different T2
nParam3 = numel(dictionaryParams(1,dictionaryParams(3,:)>0)); % number of different deviations in flip angle

nP1 = sum(dictionaryParams(1,:)>dictionaryParams(2,:));
nP2 = sum(dictionaryParams(2,:)<dictionaryParams(1,:));
nP3 = nParam3;

% Make a vector by repeating the first T1 value for as many different T2
% values there are, then the same number of the second T1 value, etc.
paramIndex1 = repmat(dictionaryParams(1,dictionaryParams(1,:)>0),nParam2, nParam3);
paramIndex1 = reshape(paramIndex1,[1,size(paramIndex1,1)*size(paramIndex1,2)*size(paramIndex1,3)]);
% Make a vector of repeats of the whole list of T1 values. Number of
% repeats equals number of different T1 values.
paramIndex2 = repmat(dictionaryParams(2,dictionaryParams(2,:)>0),[1,nParam1*nParam3]);


% normalise the dictionary
for n = 1: size(dictionary,1)
    sd(n,:) = dictionary(n,:)/dictionary(n,1);
    sd(n,:) = sd(n,:)/norm(squeeze(sd(n,:)));
end

%normalise the data
for n = 1:size(data,1)
    FPdata(n,:) = data(n,:)/data(n,1);
    FPdata(n,:) = FPdata(n,:)/norm(squeeze(FPdata(n,:)));
end


% pre-state this lines to avoid running out of memory
dictionary = reshape(dictionary,nParam1,nParam2,nParam3,nTimeCoursePts);

%% Run the matching for all voxels within the data image
tic % start timing how long the matching takes to run
disp('calculating similarity...')
for data_i = 1:(dataDim1*dataDim2) %for all pixels
    
    
    disp(['calculating similarity... fraction complete:',num2str(data_i),'/',num2str((dataDim1*dataDim2))])
    %compare the TC with each dictionary entry
    if sum(squeeze(FPdata(data_i,:))) > 0
        
        %tic
        
        timecourse = repmat(squeeze(FPdata(data_i,:)),size(sd,1),1);
        
        tempTCsimilarity = dot(timecourse',sd')';
        
        %find highest similarity score for the TC of interest
        [maxSimilarity(data_i)]  = max(tempTCsimilarity(:));
        
        %help us check if the highest similarity score found applies for
        %multiple different dictionary entries. Checks for each
  %      degeneracyMap(data_i) = sum(tempTCsimilarity(:) == max(tempTCsimilarity(:)));

        %Find the parameters associated with the best match
        tempTCsimilarity =  reshape(tempTCsimilarity',nParam1,nParam2,nParam3);
        [bestT2ind,bestT1ind,bestFAdevInd] = ind2sub(size(tempTCsimilarity),find(tempTCsimilarity == maxSimilarity(data_i)));
        
        %Assign the parameters  to the timecourse of interest
        matchedT1(data_i,1) = (dictionaryParams(1, bestT1ind));
        matchedT2(data_i,1) = (dictionaryParams(2, bestT2ind));
        matchedFAdev(data_i,1) = (dictionaryParams(3, bestFAdevInd));
        bestMatch(data_i, :) = squeeze(dictionary(bestT2ind, bestT1ind, bestFAdevInd, 1:nTimeCoursePts));
        M0fit = fit(squeeze(bestMatch(data_i, :))', data(data_i,:)',M0model,'Upper',6000,'Lower',0,'StartPoint',500 );
        M0fit_grad(data_i,1) = M0fit.a;
        
        %toc
    end
end
 
%% reshape the assigned parameter maps into 2D form
maxSimilarity = reshape(maxSimilarity,[dataDim1,dataDim2]); % Similarity score for the best match
matchedT1 = reshape(matchedT1, [dataDim1, dataDim2]);
matchedT2 = reshape(matchedT2, [dataDim1, dataDim2]);
matchedFAdev = reshape(matchedFAdev, [dataDim1, dataDim2]);
M0fit_grad = reshape(M0fit_grad, [dataDim1,dataDim2]);
bestMatch = reshape(bestMatch, [dataDim1, dataDim2,nTimeCoursePts]);

%time taken to find the best matches
matchTime = toc;

disp('Calculating similarity... Complete.');

=======
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

save([savingdir,'/MAT-files/matches/similarity/',outputPhantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'maxSimilarity.mat'],'maxSimilarity')
save([savingdir,'/MAT-files/matches/T1/',outputPhantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'matchedT1.mat'],'matchedT1')
save([savingdir,'/MAT-files/matches/T2/',outputPhantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'matchedT2.mat'],'matchedT2')
save([savingdir,'/MAT-files/matches/B1/',outputPhantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'matchedFAdev.mat'],'matchedFAdev')
save([savingdir,'/MAT-files/matches/BestMatch/',outputPhantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'bestMatch.mat'],'bestMatch')
save([savingdir,'/MAT-files/matches/MatchingTimes/',outputPhantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'compileDictionaryElapsedTime.mat'],'matchTime')
save([savingdir,'/MAT-files/matches/M0/',outputPhantomName,'offsetlist',num2str(offsetListNum),'paramList',num2str(paramRangeList),'M0fit_grad.mat'],'M0fit_grad')
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end