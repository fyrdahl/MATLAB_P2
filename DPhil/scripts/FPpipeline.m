%% Magnetic Resonance Fingerprinting Pipeline
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015

% import data to workspace
loaddata

%% process images
[bgMean,bgStd,SNR] = calc_Noise(FPimages,1,[40 2 7 7], [24 18 7 7],1);

% find signals at sample pixels
compartmentCenters = set_Compartment_Centers(phantomName);

% plot positions of sample pixels for TE and TR images
plotNumCompartments = 6;

% fit curves to calculate T1 and T2
[T1map] = fit_evolution_curves(phantomName, 'fullPhantom', 'T1', 'noPlot', TIimages, TI(2:end));

%% Simulate magnetisation evolution
%check bloch simulation by using properties of the phantom
load('fingerprintLists.mat')

% sphereD170 properties measured by preceeding fits
T1 = 282.3;
T2 = 214.8;

for TRind = 1:10
    TRs = 10:100;
    freqOffset = 0;
    nSlices = 2;
    nRepeats = 2;
    nTimeCoursePts = nRepeats*size(fingerprintLists,1);
    
    for i = 1:2
        T2s = [100,260];
        [M, Mxy(i,:),flipAngles, imageTimes, t0s] = SimSE_Bernstein(T1, T2s(i), fingerprintLists(:,:,offsetListNum), freqOffset, nSlices, nTimeCoursePts);
    end
    
    sigDiff(TRs) = sum(Mxy(1,:) - Mxy(1,:));
end

%% Create dictionary
    [signalDictionary, sdelT] = compileDictionary(fingerprintLists, offsetListNum, 48, 0, nSlices, phantomName, savingdir);

%% Check similarity scores and use dictionary to assign values of T1, T2 etc.
[similarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, match_time] = calcSimilarity(FPimages,nTimeCoursePts, dictionaryParams, paramRangeList, savingdir, phantomName, offsetListNum);
