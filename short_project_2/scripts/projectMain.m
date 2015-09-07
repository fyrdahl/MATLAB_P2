%% ONBI short project 2 - Magnetic Resonance Fingerprinting
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015
%% 1. Initialise
clear all
close all
addpath(genpath('/Applications/fsl/'))
addpath(genpath('/usr/local/fsl/bin'))
addpath(genpath('/Users/jallen/Documents/MATLAB/short_project_2'))

%% 2. Read in the images
%which phantom is the data from? ('sphereD170' or 'Jack'?)
phantomName = 'Jack';

% Choose the offset list to use. List 2 is the original 'random' list of
% offsets. Lists 3:8 are variations on List2, as described below.
%
% 3: 1st TR offset = 10000
% 4: flip Angles = 90 and 180
% 5: TE offset = 20
% 6: TR offset = 1500, TE offset = 20, FA1 = 90
% 7: TR offset = 1500, TE offset = 20
% 8: TR offset = 15000
offsetListNum = 3;

[TEImageInfo, TIImageInfo, FPImageInfo, TEimages, TIimages, FPimages, TE, TI] = readData(phantomName, offsetListNum);

%% 3. Select a ROI and a sample of the background, to calculate SNR
[SNR signal background] = calcSNR(TEimages,TE,'showFigure');

%% 4. find signals at sample pixels
run('compartmentSignals.m')

%% 5. plot positions of sample pixels for TE and TR images
plotNumCompartments = 6;
run('plotSamplePixels_TE_TR.m')
%%
run('visualiseImages.m')

%% fit curves to calculate T1 and T2

[compartmentT1s, compartmentT2s, T2curves, T1curves, fittedCurve, goodness, output] = fitEvolutionCurves(TEimages, TIimages, TE(2:end)', TI(2:end), 'compartments', compartmentCenters);

%% 6. read in the list of timing offsets used for acquisition
run('readFingerprintOffsetList.m')
fingerprintLists(:,:,offsetListNum)

%% 7. simulate magnetisation evolution
%check bloch simulation by using properties of the phantom

% sphereD170 properties
T1 = 282.5;
T2 = 214.1;

freqOffset = 0;
nSlices = 2;
[M, Mxy,flipAngles, t0s] = SimBloch(T1, T2, fingerprintLists(:,:,offsetListNum), 'showPlot', freqOffset, nSlices);

%% 8. plot positions of sample pixels for fingerprint images
compartmentCentersList = 3;
run('plotSamplePixels.m')

% plot comparison of simulation with sampled pixels
run('plotSim.m')

%% 9. create dictionary
offsetListNums = offsetListNum
run('dictionary.m')

%% 10. check similarity and use dictionary to measure T1 and T2
run('matching.m')