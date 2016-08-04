function [sd,sdelT] = compileDictionary(fingerprintLists, offsetListNum, nTimeCoursePts, freqOffset, nSlices, phantomName, savingdir)

%% DICTIONARY

%% Import fingerprint timings list and dictionary parameters
load /Users/jallen/Documents/short_project_2/MAT-files/fingerprintLists.mat
[dictionaryParams, paramList] = setDictionaryParams(phantomName,3);

%%
sd = zeros(sum(dictionaryParams(1,:)>0), sum(dictionaryParams(2,:)>0), sum(dictionaryParams(3,:)>0) , nTimeCoursePts, max(offsetListNum));

originalFA1s = fingerprintLists(:,3,offsetListNum);
originalFA2s = fingerprintLists(:,4,offsetListNum);


tic
%vary T1
for i = 1:sum(dictionaryParams(1,:)>0)

    T1 = dictionaryParams(1,i);
    % vary T2
    for j = 1:sum(dictionaryParams(2,:)>0)
        T2 = dictionaryParams(2,j);
        
        % vary intended flip angle, to allow for excitation field inhomogeneties 
        for k = 1:sum(dictionaryParams(3,:)>0)
            % apply B1 variation range
            fingerprintLists(:,3,offsetListNum) = originalFA1s*(dictionaryParams(3,k));
            fingerprintLists(:,4,offsetListNum) = originalFA2s*(dictionaryParams(3,k));

            if T1>=T2 % T1 must be larger than T2
            [~, sd(i,j,k,:,offsetListNum), ~, ~] = SimSE_Bernstein(T1, T2, fingerprintLists(:,:,offsetListNum), freqOffset, nSlices, nTimeCoursePts);
            end
            %% add noise to the simulated signals
            %             SNR = 0.655*squeeze(signalDictionary(i,j,k,:))/std(background(:));
            %             for tPt = 1:size(signalDictionary,4)
            %                 signalDictionary(i,j,k,tPt) =  awgn(signalDictionary(i,j,k,tPt),SNR(tPt));
            %             end
            
        end
    end
    disp(['compiling dictionary for offset list ', num2str(offsetListNum),': ',num2str( (i/(sum(dictionaryParams(1,:)~=0)) )*100) , ' percent complete'])

end
toc
sdelT = toc;
%% Save the dictionary
save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'dictionary.mat'],'sd')
save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'signalDictionaryTime.mat'],'sdelT')

end


