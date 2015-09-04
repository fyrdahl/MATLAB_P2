%% ONBI short project 2 - Magnetic Resonance Fingerprinting
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015

clear all
close all
addpath(genpath('/Applications/fsl/'))
addpath(genpath('/usr/local/fsl/bin'))
addpath(genpath('/Users/jallen/Documents/MATLAB/short_project_2'))

%% Read in the images
%which phantom is the data from? ('sphereD170' or 'Jack'?)
phantomName = 'sphereD170';

% Choose the offset list to use. List 2 is the original 'random' list of
% offsets. Lists 3:8 are variations on List2, as described below.
%
% 3: 1st TR offset = 10000
% 4: flip Angles = 90 and 180
% 5: TE offset = 20
% 6: TR offset = 1500, TE offset = 20, FA1 = 90
% 7: TR offset = 1500, TE offset = 20
% 8: TR offset = 15000
offsetListNum = 7;

[TEImageInfo, TIImageInfo, FPImageInfo, TEimages, TIimages, FPimages, TE, TI] = readData(phantomName, offsetListNum);

%% select a ROI and a sample of the background, to calculate SNR
[SNR signal background] = calcSNR(TEimages,TE,'showFigure');

%%
run('compartmentSignals.m')

%%
run('visualiseImages.m')

%% fit curves to calculate T1 and T2

[compartmentT1s, compartmentT2s, T2curves, T1curves, fittedCurve, goodness, output] = fitEvolutionCurves(TEimages, TIimages, TE(2:end)', TI(2:end), 'compartments', compartmentCenters);

%% read in the list of timing offsets used for acquisition
run('readFingerprintOffsetList.m')
fingerprintLists(:,:,offsetListNum)

%% simulate magnetisation evolution
%check bloch simulation by using properties of the phantom

% sphereD170 properties
T1 = 282.5;
T2 = 214.1;

freqOffset = 0;
nSlices = 2;
[M, Mxy,flipAngles, t0s] = SimBloch(T1, T2, fingerprintLists(:,:,offsetListNum), 'showPlot', freqOffset, nSlices);

%
plotNumCompartments = 6;
%% Show the positions of the sampled compartment centers, with respect to the images
figure;
imagesc(TEimages(:,:,TE(1)))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1 :plotNumCompartments
    plot(compartmentCenters(i,1,1),compartmentCenters(i,2,1),'*')
    text(compartmentCenters(i,1,1),compartmentCenters(i,2,1), compartmentLabels(i) )
end

figure;
imagesc(TIimages(:,:,TI(1)))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1:plotNumCompartments
    plot(compartmentCenters(i,1,2),compartmentCenters(i,2,2),'*')
    text(compartmentCenters(i,1,2),compartmentCenters(i,2,2), compartmentLabels(i) )
end

figure;
imagesc(FPimages(:,:,1))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1:plotNumCompartments
    plot(compartmentCenters(i,1,3),compartmentCenters(i,2,3),'*')
    text(compartmentCenters(i,1,3),compartmentCenters(i,2,3), compartmentLabels(i) )
end
%% plot positions of sample pixels
compartmentCentersList = 1;

figure;
imagesc(FPimages(:,:,1))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1:plotNumCompartments
    plot(compartmentCenters(i,1,compartmentCentersList),compartmentCenters(i,2,compartmentCentersList),'*')
    text(compartmentCenters(i,1,compartmentCentersList),compartmentCenters(i,2,compartmentCentersList), compartmentLabels(i) )
end
%% plot comparison of simulation with sampled pixels
run('plotSim.m')

%%

% normalise dictionary entries to have the same sum squared magnitude
% select one dictionary entry for each pixel, using the complex data for simulation and pixel
% calculate proton density (M0) as the scaling factor between the measured
% signal and the simulated (see Ma2013).

%% create dictionary
clear dictionaryParams
dictionaryParams(1,:) = 200:10:300 ; % T1
dictionaryParams(2,:) = 200:10:300 ; % T2
dictionaryParams(3,:) = 0.8:0.04:1.2 ; % B1 fraction

nTimeCoursePts = size(data , 4)/2;

signalDictionary = zeros(size(dictionaryParams(1,:),2), size(dictionaryParams(2,:),2), size(dictionaryParams(3,:),2), nTimeCoursePts);

for offsetListNum = 2:8
[signalDictionary(:,:,:,:,offsetListNum)] = compileDictionary(fingerprintLists, offsetListNum, dictionaryParams, nTimeCoursePts, freqOffset, nSlices, background);
% !!! must normalise dictionary entries to have the same sum squared
% magnitude (use sumsqr() )
end

%% check similarity and use dictionary to measure T1 and T2
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

[similarity, matchedT1, matchedT2, matchedFAdevInd, bestMatch] = calcSimilarity(data, signalDictionary(:,:,:,:,offsetListNum), sliceNumber, dictionaryParams);
%%
%%visualise spread of matched T1s and T2s
figure; hist(squeeze(matchedT1(offsetListNum,:)))
figure;
hist(squeeze(matchedT2(offsetListNum,:)))
%%
figure; plot(squeeze(data(1,1,1,:)), '-*')
hold on
plot(squeeze(bestMatch(1,1,:))*data(1,1,1,1),'--.')
