function [dictionary,sdelT] = compile_SE_dictionary_EPG(offsetList, nPts, nSlices, df, dictParams,TRmin, TEmin)
%% Jack Allen
%
%% DICTIONARY



%%
tic

originalFA1s = offsetList(:,3);
originalFA2s = offsetList(:,4);
%
paramIndex1 = repmat(dictParams(1,dictParams(1,:)>0),sum(dictParams(2,:)>0), sum(dictParams(3,:)>0));
paramIndex1 = reshape(paramIndex1,1,size(paramIndex1,1)*size(paramIndex1,2)*size(paramIndex1,3));
paramIndex2 = repmat(dictParams(2,dictParams(2,:)>0),1,sum(dictParams(1,:)>0)*sum(dictParams(3,:)>0));
paramIndex3 = repmat(dictParams(3,dictParams(3,:)>0),sum(dictParams(1,:)>0)*sum(dictParams(2,:)>0),1);
paramIndex3 = reshape(paramIndex3,1,size(paramIndex3,1)*size(paramIndex3,2));

dictionary = zeros(size(paramIndex1,2), nPts);

%vary T1
disp(['building signal dictionary...(Number of entries:',num2str(numel(find(paramIndex1>=paramIndex2))),').'])


for nEntry = find(paramIndex1>paramIndex2)
    % vary intended flip angle, to allow for excitation field inhomogeneties
    offsetList(:,3) = paramIndex3(nEntry)*originalFA1s; %FA1 variation
    offsetList(:,4) = paramIndex3(nEntry)*originalFA2s; %FA2 variation
    
    [dictionary(nEntry,:)] = sim_SE_EPG(paramIndex1(nEntry), paramIndex2(nEntry), TRmin, TEmin, offsetList,(nPts/24), nSlices, df, 'no plot');  
end

disp('building signal dictionary: complete')
toc
sdelT = toc;

end


