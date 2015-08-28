function [signalDictionary] = compileDictionary(fingerprintLists, offsetListNum, dictionaryParams, nTimeCoursePts, freqOffset, nSlices)

%% DICTIONARY

nTimeCoursePts = 24;

freqOffset = 0;

nSlices = 2;

signalDictionary = zeros(size(dictionaryParams(1,:),2), size(dictionaryParams(2,:),2), nTimeCoursePts);

tic 
for i = 1:numel(dictionaryParams(1,:))
    
    T1 = dictionaryParams(1,i);
      
    for j = 1:numel(dictionaryParams(2,:))
        
        T2 = dictionaryParams(2,j);
        
        [~, ~, ~,  signalDictionary(i,j,:), ~, ~] = SimBloch(T1, T2, fingerprintLists(:,:,offsetListNum), 'dontPlot', freqOffset, nSlices);
         
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
