function [matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime, maxSimilarity, degeneracyMap] = calc_similarity(data, dictionary, dictionaryParams)
% [maxSimilarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime] = calc_similarity(data, nTimeCoursePts, dictionaryParams, paramRangeList,savingdir, inputPhantomName,outputPhantomName,offsetListNum)
% Description: A function to calculate the similarity of a time course with
% each entry in a dictionary of time courses and find the best match.
%
% Author: Jack Allen.
% Institution: University of Oxford.
% Contact: jack.allen@jesus.ox.ac.uk
%
%
%% Similarity measurements for magnetic resonance fingerprinting (MRF) matching
disp('calculating similarity: started')

%mask to exclude time course outside of ROI
%loadMask.Mask = load([savingdir,'/MAT-files/mask.mat']);
%mask = reshape(loadMask.mask,[64*64, 1]);

nTimeCoursePts = size(data,3);

dataDim1 = size(data,1);
dataDim2 = size(data,2);
%time course of interest
data = reshape(data,[size(data,1)*size(data,2), nTimeCoursePts]);

%% Initialise maps
maxSimilarity = zeros(size(data,1),1);
degeneracyMap = zeros(size(data,1),1);
bestMatch = zeros(size(data,1),nTimeCoursePts);
matchedT1 = zeros(size(data,1),1);
matchedT2 = zeros(size(data,1),1);
matchedFAdev = zeros(size(data,1),1);
M0fit_grad = zeros(size(data,1),1);

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

end