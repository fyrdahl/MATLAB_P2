%% define data
sliceNumber = 1 % slice to be analysed
clear data

% data = FPimages(compartmentCenters(1,1),compartmentCenters(1,2),:,:);

for r = 1:size(compartmentCenters,1)
    data(r,1,sliceNumber,:) = FPimages(compartmentCenters(r,1),compartmentCenters(r,2),sliceNumber,:);
end

% for r = 1 : size(FPimages,1)
%     for c = 1 : size(FPimages,2)
% data(r,c,sliceNumber,:) = FPimages(r,c,sliceNumber,:);
%     end
% end

%% create dictionary

switch phantomName
    case 'sphereD170'
        
        disp('Phantom: sphereD170')
        clear dictionaryParams
        dictionaryParams(1,:) = 200:10:300 ; % T1
        dictionaryParams(2,:) = 200:10:300 ; % T2
        dictionaryParams(3,:) = 0.7:0.06:1.3 ; % B1 fraction
    case 'Jack'
        disp('Phantom: Jack')
        clear dictionaryParams
        dictionaryParams(1,:) = [100:50:200,400:50:500,800:50:900,2000:50:2100,2500:50:3000] ; % T1
        dictionaryParams(2,:) = [20:(250-20)/(numel(dictionaryParams(1,:))-1):250] ; % T2
        dictionaryParams(3,:) = [0.7:((1.3 - 0.7)/(numel(dictionaryParams(1,:))-1)):1.3] ; % B1 fraction
end
nTimeCoursePts = size(data , 4)/2;

signalDictionary = zeros(size(dictionaryParams(1,:),2), size(dictionaryParams(2,:),2), size(dictionaryParams(3,:),2), nTimeCoursePts);


[signalDictionary] = compileDictionary(fingerprintLists, offsetListNums, dictionaryParams, nTimeCoursePts, freqOffset, nSlices, background);
% !!! must normalise dictionary entries to have the same sum squared
% magnitude (use sumsqr() )



% normalise dictionary entries to have the same sum squared magnitude
% select one dictionary entry for each pixel, using the complex data for simulation and pixel
% calculate proton density (M0) as the scaling factor between the measured
% signal and the simulated (see Ma2013).

