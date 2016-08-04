% Magnetic Resonance Fingerprinting Pipeline
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015

% import data to workspace
loadFPdata

%% process images
[bgMean,bgStd,SNR] = calc_noise(FPimages,1,[40 2 7 7], [24 18 7 7],0);

% fit curves to calculate T1 and T2
%[T1map] = fit_evolution_curves(phantomName, 'fullPhantom', 'T1', 'noPlot', TIimages, TI(2:end));

%% Simulate magnetisation evolution
%check bloch simulation by using properties of the phantom
load('fingerprintLists.mat')

% sphereD170 properties measured by preceeding fits
T1 = 282.3;
T2 = 214.8;

offsetListNum = 3;
freqOffset = 0;
nSlices = 2; % How many slices?
nRepeats = 2; % How many cycles through the offset list?
nTimeCoursePts = nRepeats*size(fingerprintLists,1); % How many timepoints

[M, Mxy,flipAngles, imageTimes, t0s] = SimSE_Bernstein(T1, T2, fingerprintLists(:,:,offsetListNum), freqOffset, nSlices, nTimeCoursePts);

%% Create dictionary
paramList = 2; % specify the list of dictionary parameters
dictionaryParams = set_dictionary_params('Jack1',4);


[dictionary,sdelT] = compile_SE_dictionary(offsetList, nTimeCoursePts, freqOffset, dictionaryParams);

% Check similarity scores and use dictionary to assign values of T1, T2 etc.
[similarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, match_time] = calc_similarity(FPimages,nTimeCoursePts, paramList, savingdir, phantomName, offsetListNum);
