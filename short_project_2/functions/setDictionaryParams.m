function [dictionaryParams] = setDictionaryParams(phantomName)

switch phantomName
    case 'sphereD170'
        
        disp('Phantom: sphereD170')
    
        T1s = 200:10:300
        T2s = 200:10:300
        FAdevs = 0.7:0.05:1.3
        dictionaryParams(1,1:numel(T1s)) = T1s; % T1
        dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
        
    case 'Jack'
        
        disp('Phantom: Jack')
     
        T1s = [30:20:260, 3050:20:3150]
        T2s = [10:10:120, 1770:10:1820]
        FAdevs = 0.8:0.05:1.2
        dictionaryParams(1,1:numel(T1s)) = T1s;
        dictionaryParams(2,1:numel(T2s)) = T2s;
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs;
        
end
end
