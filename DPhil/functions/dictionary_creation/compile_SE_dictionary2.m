function [sd,sdelT] = compile_SE_Dictionary2(fingerprintLists, offsetListNum, nTimeCoursePts, freqOffset, nSlices, phantomName, savingdir)

%% DICTIONARY

%% Import fingerprint timings list and dictionary parameters
load /Users/jallen/Documents/short_project_2/MAT-files/fingerprintLists.mat
[dictionaryParams, paramList] = setDictionaryParams(phantomName,3);

%%

sd = zeros(sum(dictionaryParams(1,:)>0)*sum(dictionaryParams(2,:)>0)*sum(dictionaryParams(3,:)>0) , nTimeCoursePts);
originalFA1s = fingerprintLists(:,3,offsetListNum);
originalFA2s = fingerprintLists(:,4,offsetListNum);
paramIndex1 = repmat(dictionaryParams(1,dictionaryParams(1,:)>0),sum(dictionaryParams(2,:)>0), sum(dictionaryParams(3,:)>0));
paramIndex1 = reshape(paramIndex1,1,size(paramIndex1,1)*size(paramIndex1,2)*size(paramIndex1,3));
paramIndex2 = repmat(dictionaryParams(2,dictionaryParams(2,:)>0),1,sum(dictionaryParams(1,:)>0)*sum(dictionaryParams(3,:)>0));
paramIndex3 = repmat(dictionaryParams(3,dictionaryParams(3,:)>0),sum(dictionaryParams(1,:)>0)*sum(dictionaryParams(2,:)>0),1);
paramIndex3 = reshape(paramIndex3,1,size(paramIndex3,1)*size(paramIndex3,2));

%remove unphysical combinations of T1 and T2
tempparamIndex1 = paramIndex1;
tempparamIndex2 = paramIndex2;
paramIndex1(tempparamIndex1<tempparamIndex2) = [];
paramIndex2(tempparamIndex1<tempparamIndex2) = [];
paramIndex3(tempparamIndex1<tempparamIndex2) = [];
nDictionaryEntries = size(paramIndex3,2);

tic
%vary T1
disp('building signal dictionary...')
for i = 1:nDictionaryEntries
    
    % vary intended flip angle, to allow for excitation field inhomogeneties
    fingerprintLists(:,3,offsetListNum) = originalFA1s*paramIndex3(i);
    fingerprintLists(:,4,offsetListNum) = originalFA2s*paramIndex3(i);
    [~, sd(i,:), ~, ~] = SimSE_Bernstein(paramIndex1(i), paramIndex2(i), fingerprintLists(:,:,offsetListNum), freqOffset, nSlices, nTimeCoursePts);
    
    %% add noise to the simulated signals
    %             SNR = 0.655*squeeze(signalDictionary(i,j,k,:))/std(background(:));
    %             for tPt = 1:size(signalDictionary,4)
    %                 signalDictionary(i,j,k,tPt) =  awgn(signalDictionary(i,j,k,tPt),SNR(tPt));
    %             end
end
disp('building signal dictionary: complete')
toc
sdelT = toc;
%% Save the dictionary
save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'dictionary.mat'],'sd')
save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'signalDictionaryTime.mat'],'sdelT')

end


