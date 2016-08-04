<<<<<<< HEAD
function [dictionary,sdelT] = compile_SE_dictionary(offsetList, nTimeCoursePts, nSlices, freqOffset, dictionaryParams)
%% Jack Allen
%
%% DICTIONARY



%%
tic

originalFA1s = offsetList(:,3);
originalFA2s = offsetList(:,4);
%
paramIndex1 = repmat(dictionaryParams(1,dictionaryParams(1,:)>0),sum(dictionaryParams(2,:)>0), sum(dictionaryParams(3,:)>0));
paramIndex1 = reshape(paramIndex1,1,size(paramIndex1,1)*size(paramIndex1,2)*size(paramIndex1,3));
paramIndex2 = repmat(dictionaryParams(2,dictionaryParams(2,:)>0),1,sum(dictionaryParams(1,:)>0)*sum(dictionaryParams(3,:)>0));
paramIndex3 = repmat(dictionaryParams(3,dictionaryParams(3,:)>0),sum(dictionaryParams(1,:)>0)*sum(dictionaryParams(2,:)>0),1);
paramIndex3 = reshape(paramIndex3,1,size(paramIndex3,1)*size(paramIndex3,2));

dictionary = zeros(size(paramIndex1,2), nTimeCoursePts);

%vary T1
disp(['building signal dictionary...(Number of entries:',num2str(numel(find(paramIndex1>=paramIndex2))),').'])

TRmin = 130;
TEmin = 32;

for nEntry = find(paramIndex1>paramIndex2)
    % vary intended flip angle, to allow for excitation field inhomogeneties
    offsetList(:,3) = paramIndex3(nEntry)*originalFA1s; %FA1 variation
    offsetList(:,4) = paramIndex3(nEntry)*originalFA2s; %FA2 variation
    
    %[Msig, M, TRindex] = sim_SE2(T1, T2, fingerprintOffsetList,nRepeats,df, plotFlag)
    [dictionary(nEntry,:)] = sim_SE3(paramIndex1(nEntry), paramIndex2(nEntry), TRmin, TEmin, offsetList,(nTimeCoursePts/24), nSlices, freqOffset, 'no plot');
    %[~, dictionary(nEntry,:)] = sim_SE_bernstein(paramIndex1(nEntry), paramIndex2(nEntry), offsetList,freqOffset, 2, nTimeCoursePts,'noPlot');
    %% add noise to the simulated signals
    %             SNR = 0.655*squeeze(signalDictionary(i,j,k,:))/std(background(:));
    %             for tPt = 1:size(signalDictionary,4)
    %                 signalDictionary(i,j,k,tPt) =  awgn(signalDictionary(i,j,k,tPt),SNR(tPt));
    %             end
    
end

disp('building signal dictionary: complete')
toc
sdelT = toc;
=======
function [sd,sdelT] = compile_SE_Dictionary(fingerprintLists, offsetListNum, nTimeCoursePts, freqOffset, nSlices, phantomName, savingdir)

%% DICTIONARY

%% Import fingerprint timings list and dictionary parameters
load /Users/jallen/Documents/short_project_2/MAT-files/fingerprintLists.mat
[dictionaryParams, paramList] = setDictionaryParams(phantomName,3);

%%
sd = zeros(sum(dictionaryParams(1,:)>0), sum(dictionaryParams(2,:)>0), sum(dictionaryParams(3,:)>0) , nTimeCoursePts, 9);

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
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1

end


