%% define data
sliceNumber = 1 % slice to be analysed
clear data

% data = FPimages(compartmentCenters(1,1),compartmentCenters(1,2),:,:);

% for r = 1:size(compartmentCenters,1)
%     data(r,1,sliceNumber,:) = FPimages(compartmentCenters(r,1),compartmentCenters(r,2),sliceNumber,:,offsetListNum);
% end

for r = 1 : size(FPimages,1)
    for c = 1 : size(FPimages,2)
data(r,c,sliceNumber,:) = FPimages(r,c,sliceNumber,:);
    end
end

%% create dictionary

switch phantomName
    case 'sphereD170'
        
        disp('Phantom: sphereD170')
        clear dictionaryParams
        T1s = 200:10:300;
        T2s = 200:10:300;
        FAdevs = 0.7:0.05:1.3 ;
        
        dictionaryParams(1,1:numel(T1s)) = T1s;
        dictionaryParams(2,1:numel(T2s)) = T2s;
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs;
        
    case 'Jack'
        
        disp('Phantom: Jack')
        clear dictionaryParams
        T1s = [30:20:260, 3050:20:3150];
        T2s = [10:10:120, 1770:10:1820];
        FAdevs = 0.7:0.05:1.3;
        
        dictionaryParams(1,1:numel(T1s)) = T1s;
        dictionaryParams(2,1:numel(T2s)) = T2s;
        dictionaryParams(3,1:numel(FAdevs)) = FAdevs;
        
end
nTimeCoursePts = 24;


[signalDictionary] = compileDictionary(fingerprintLists, offsetListNums, dictionaryParams, nTimeCoursePts, freqOffset, nSlices, background);
% !!! must normalise dictionary entries to have the same sum squared
% magnitude (use sumsqr() )

% normalise dictionary entries to have the same sum squared magnitude
% select one dictionary entry for each pixel, using the complex data for simulation and pixel
% calculate proton density (M0) as the scaling factor between the measured
% signal and the simulated (see Ma2013).

