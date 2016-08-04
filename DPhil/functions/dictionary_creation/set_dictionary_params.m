<<<<<<< HEAD
function [dictionaryParams, paramList] = set_dictionary_params(phantomName,varargin)
=======
function [dictionaryParams, paramList] = setDictionaryParams(phantomName,paramList)
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1

switch phantomName
    case 'sphereD170'
        
        disp('Phantom: sphereD170')
<<<<<<< HEAD
        switch varargin{1}
=======
        switch paramList
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
            case 1
                paramList = 1;
                T1s = 200:10:300;
                T2s = 200:10:300;
                FAdevs = 0.7:0.05:1.3;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
<<<<<<< HEAD
            case 2
                paramList = 2;
                T1s = 275:1:290;
                T2s = 210:1:220;
                FAdevs = 0.7:0.05:1.3;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
=======
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
            case 3
                paramList = 3;
                T1s = 200:10:300;
                T2s = 200:10:300;
                FAdevs = 0.7:0.05:1.3;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
<<<<<<< HEAD
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
            case 4
                paramList = 4;
                T1s = 200:5:300;
                T2s = 200:5:300;
                FAdevs = 1;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
            case 5
                paramList = 5;
                T1s = 200:10:300;
                T2s = 200:10:300;
                FAdevs = 0.7:0.005:1.3;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
                
            case 6
                paramList = 6;
                T1s = 200:10:300;
                T2s = 200:10:300;
                FAdevs = 0.7:0.01:1.3;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
                
            case 7
                paramList = 7;
                T1s = 200:10:300;
                T2s = 200:10:300;
                FAdevs = 1;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
            case 8
                paramList = 8;
                T1s = 200:5:300;
                T2s = 200:5:300;
                FAdevs =  0.7:0.01:1.3;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
                
            case 9
                paramList = 9;
                T1s = 200:5:300;
                T2s = 200:5:300;
                FAdevs =  0.5:0.01:1.5;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction

            case 10
                paramList = 10;
                T1s = 200:1:300;
                T2s = 200:1:300;
                FAdevs =  0.7:0.01:1.7;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
                
                 case 11
                paramList = 11;
                T1s = 200:1:300;
                T2s = 200:1:300;
                FAdevs =  1;
                dictionaryParams(1,1:numel(T1s)) = T1s; % T1
                dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
                
        end
        
    case 'Jack1'
        disp('Phantom: Jack1')
        switch varargin{1}
=======
                dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction       
        end
            
    case 'Jack'
        disp('Phantom: Jack')
        switch paramList
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
            case 1
                paramList = 1;
                T1s = [30:10:300, 3050:10:3200];
                T2s = [10:10:120, 1900:10:2100];
                FAdevs = 0.7:0.05:1.3;
            case 2
                paramList = 2;
                T1s = [30:20:260, 3050:20:3150];
                T2s = [10:10:120, 1770:10:1820];
                FAdevs = 0.7:0.05:1.3;
            case 3
                paramList = 3;
                T1s = [30:20:260, 3050:20:3150];
                T2s = [10:10:120, 1770:10:1820];
                FAdevs = 0.7:0.01:1.3;
<<<<<<< HEAD
            case 4
                paramList = 4;
                T1s = [20:10:300, 3000:10:3200];
                T2s = [10:10:150, 1750:10:2000];
                FAdevs = 0.7:0.01:1.3;
=======
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
        end
        
        dictionaryParams(1,1:numel(T1s)) = T1s;
        dictionaryParams(2,1:numel(T2s)) = T2s;
<<<<<<< HEAD
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs;
        
    case 'brainTemplate'
        
        paramList = 'brainTemplate';
        T1s = [750:5:850,1150:5:1250,2950:5:3050]; % around 800ms, 1200ms, 3000ms
        T2s = [40:5:90,450:5:550]; % around 60ms, 75ms, 500ms
        FAdevs = 0.7:0.005:1.3;
        dictionaryParams(1,1:numel(T1s)) = T1s; % T1
        dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction
        
=======
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs;    
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end
end
