function [dictionaryParams, paramList] = setDictionaryParams(phantomName,paramList)

switch phantomName
    case 'sphereD170'
        
        disp('Phantom: sphereD170')
     switch paramList
         case 1
        paramList = 1
        T1s = 200:10:300;
        T2s = 200:10:300;
        FAdevs = 0.7:0.05:1.3;
        dictionaryParams(1,1:numel(T1s)) = T1s; % T1
        dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
    case 3
        paramList = 3
        T1s = 200:10:300;
        T2s = 200:10:300;
        FAdevs = 0.7:0.01:1.3;
        dictionaryParams(1,1:numel(T1s)) = T1s; % T1
        dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
   
     end
     
     
    case 'Jack'
        
        disp('Phantom: Jack')
     switch paramList
         case 1
             paramList = 1
        T1s = [30:10:300, 3050:10:3200]
        T2s = [10:10:120, 1900:10:2100]       
         FAdevs = 0.7:0.05:1.3
         case 2
             paramList = 2
        T1s = [30:20:260, 3050:20:3150]
        T2s = [10:10:120, 1770:10:1820]
        FAdevs = 0.7:0.05:1.3
        case 3
             paramList = 3
        T1s = [30:20:260, 3050:20:3150]
        T2s = [10:10:120, 1770:10:1820]
        FAdevs = 0.7:0.01:1.3
     end
        
        dictionaryParams(1,1:numel(T1s)) = T1s;
        dictionaryParams(2,1:numel(T2s)) = T2s;
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs;
        
end
end
