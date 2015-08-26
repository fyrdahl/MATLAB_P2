%% DICTIONARY
T1dictionaryList = [50;100;150;200;250;300]
T2dictionaryList = [50;100;150;200;250;300]
for T1 = T1dictionaryList'
    T1
    for T2 = T2dictionaryList'
        T2
        [Mtransverse(find(T1dictionaryList == T1),find(T2dictionaryList == T2)) M] = SimBloch(T1, T2, fingerprintList2, 'showPlot');
    end
end
