function [dictionary,paramIndices, sdelT] = compile_SE_dictionary_Bernstein(offsets, minTR, minTE, nRuns, df, dictionaryParams)
%% Jack Allen
%
%% DICTIONARY



%%
tic

originalFA1s = offsets(:,3);
originalFA2s = offsets(:,4);
%
paramIndex1 = repmat(dictionaryParams(1,dictionaryParams(1,:)>0),sum(dictionaryParams(2,:)>0), sum(dictionaryParams(3,:)>0));
paramIndex1 = reshape(paramIndex1,1,size(paramIndex1,1)*size(paramIndex1,2)*size(paramIndex1,3));
paramIndex2 = repmat(dictionaryParams(2,dictionaryParams(2,:)>0),1,sum(dictionaryParams(1,:)>0)*sum(dictionaryParams(3,:)>0));
paramIndex3 = repmat(dictionaryParams(3,dictionaryParams(3,:)>0),sum(dictionaryParams(1,:)>0)*sum(dictionaryParams(2,:)>0),1);
paramIndex3 = reshape(paramIndex3,1,size(paramIndex3,1)*size(paramIndex3,2));

paramIndices(1,:) = paramIndex1;
paramIndices(2,:) = paramIndex2;
paramIndices(3,:) = paramIndex3;

dictionary = zeros(size(paramIndex1,2), size(offsets,1)*nRuns);

disp(['building signal dictionary...(Number of entries:',num2str(numel(find(paramIndex1>=paramIndex2))),').'])
for nEntry = find(paramIndex1>paramIndex2)
    % vary intended flip angle, to allow for excitation field inhomogeneties
    offsets(:,3) = paramIndex3(nEntry)*originalFA1s; %FA1 variation
    offsets(:,4) = paramIndex3(nEntry)*originalFA2s; %FA2 variation

   % Msig = sim_SE_bernstein(T1, T2, TRmin, TEmin, offsets,df, nSlices, size(data,3));
   dictionary(nEntry,:) = sim_SE_bernstein(paramIndex1(nEntry), paramIndex2(nEntry), minTR, minTE, offsets,nRuns,df);
    %% add noise to the simulated signals
    %             SNR = 0.655*squeeze(signalDictionary(i,j,k,:))/std(background(:));
    %             for tPt = 1:size(signalDictionary,4)
    %                 signalDictionary(i,j,k,tPt) =  awgn(signalDictionary(i,j,k,tPt),SNR(tPt));
    %             end
    
end

disp('building signal dictionary: complete')
toc
sdelT = toc;

end


