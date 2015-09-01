function [signalDictionary] = compileDictionary(fingerprintLists, offsetListNum, dictionaryParams, nTimeCoursePts, freqOffset, nSlices)

%% DICTIONARY
originalFAs = fingerprintLists(:,3,offsetListNum);

signalDictionary = zeros(size(dictionaryParams(1,:),2), size(dictionaryParams(2,:),2), size(dictionaryParams(3,:),2) , nTimeCoursePts);

tic 
for i = 1:numel(dictionaryParams(1,:))
    
    %vary T1
    T1 = dictionaryParams(1,i);
      
    for j = 1:numel(dictionaryParams(2,:))
        
        % vary T2
        T2 = dictionaryParams(2,j);
        
        for k = 1:numel(dictionaryParams(3,:))
            
         % vary flip angle 1
        fingerprintLists(:,3,offsetListNum) = originalFAs + originalFAs.*dictionaryParams(3,k);
              
%         [~, ~, ~,  signalDictionary(i,j,k,:), ~, ~] =
        a = SimBloch(T1, T2, fingerprintLists(:,:,offsetListNum), 'dontPlot', freqOffset, nSlices);
        
        end
    end
    
    disp(['compiling dictionary for list ', num2str(offsetListNum),': ',num2str( (i/numel(dictionaryParams(1,:)))*100) , ' percent complete'])
end
toc
%%
% figure;
% hold
% plot(Mxy,'-.*')
% for i = 1:numel(dictionaryParams(1,:))
%     for j = 1:numel(dictionaryParams(2,:))
% plot(squeeze(signalDictionary(i,j,:)),'-^')
%     end
% end

end
