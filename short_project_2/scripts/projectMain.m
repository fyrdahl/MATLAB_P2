%% ONBI short project 2 - Magnetic Resonance Fingerprinting
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015
%% 1. Initialise
clear all
close all

%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/short_project_2';
workingdir = '/Users/jallen/Documents/MATLAB/short_project_2';

addpath(genpath(workingdir));
addpath(genpath('/Applications/fsl/'))
addpath(genpath('/usr/local/fsl/bin'))
addpath(genpath('/opt/fmrib/fsl/etc/matlab'))

%% 2. Read in the images
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

for offsetListNum = 2:8
[TEImageInfo, TIImageInfo, FPImageInfo(:,offsetListNum), TEimages, TIimages, FPimages(:,:,:,:,offsetListNum), TE, TI] = readData(phantomName, offsetListNum, workingdir);
save([workingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TEimages.mat'],'TEimages')
save([workingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'TIimages.mat'],'TIimages')
save([workingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'FPimages.mat'],'FPimages')
end
%% 3. Select a ROI and a sample of the background, to calculate SNR
[SNR signal background] = calcSNR(TEimages,TE,'showFigure');

%% 4. find signals at sample pixels
run('compartmentSignals.m')

%% 5. plot positions of sample pixels for TE and TR images
plotNumCompartments = 6;
sliceNumber = 1;
run('plotSamplePixels_TE_TR.m')
%%
run('visualiseImages.m')

%% fit curves to calculate T1 and T2

[compartmentT1s, compartmentT2s, T2curves, T1curves, fittedCurve, goodness, output] = fitEvolutionCurves(TEimages, TIimages, TE(2:end)', TI(2:end), 'compartments', compartmentCenters);

%% 6. read in the list of timing offsets used for acquisition
run('readFingerprintOffsetList.m')
save([workingdir,'/MAT-files/fingerprintLists.mat'], 'fingerprintLists')
%% 7. simulate magnetisation evolution
%check bloch simulation by using properties of the phantom
load('fingerprintLists.mat')
% sphereD170 properties
T1 = 282.5;
T2 = 214.1;

freqOffset = 0;
nSlices = 2;
[M, Mxy,flipAngles, t0s] = SimBloch(T1, T2, fingerprintLists(:,:,offsetListNum), 'showPlot', freqOffset, nSlices);

%% 8. check signal simulation by plotting positions of sample pixels for the fingerprinting images
compartmentCentersList = 3;
run('plotSamplePixels.m')

% plot comparison of simulation with sampled pixels
run('plotSim.m')

%% 9. create dictionary
load fingerprintLists.mat

dictionaryParams = setDictionaryParams(phantomName);

nTimeCoursePts = 24;
for offsetListNum = 2:8;
[signalDictionary] = compileDictionary(fingerprintLists, offsetListNum, dictionaryParams, nTimeCoursePts, freqOffset, nSlices, phantomName, background);
save([workingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'dictionary.mat'],'signalDictionary')
pause(1)
end

%% 10. check similarity and use dictionary to measure T1 and T2

sliceNumber = 1 ;% slice to be analysed
clear data

% data = FPimages(compartmentCenters(1,1),compartmentCenters(1,2),:,:);

% for r = 1:size(compartmentCenters,1)
% data(r,1,sliceNumber,:) = FPimages(compartmentCenters(r,1),compartmentCenters(r,2),sliceNumber,:);
% end
%%
for offsetListNum = 2:8;
clear FPimages
clear data
load([workingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'FPimages.mat'])
data = zeros(size(FPimages,1), size(FPimages,2), 1, 24);
for r = 1 : size(FPimages,1)
    for c = 1 : size(FPimages,2)
data(r,c,sliceNumber,1:24) = squeeze(FPimages(r,c,sliceNumber,1:24,offsetListNum));
    end
end

    disp(['offsetList',num2str(offsetListNum)])
clear signalDictionary
l = load([workingdir,'/MAT-files/dictionaries/list',num2str(offsetListNum),'dictionary.mat']);
signalDictionary = l.signalDictionary;
[similarity(:,:,:,:,:), matchedT1(:,:), matchedT2(:,:), matchedFAdevInd(:,:)] = calcSimilarity(data, signalDictionary(:,:,:,:,offsetListNum), sliceNumber, dictionaryParams, workingdir);
save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'matchedT1.mat'],'matchedT1')
save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'matchedT2.mat'],'matchedT2')
save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'matchedFAdevInd.mat'],'matchedFAdevInd')
pause(1)

end
%
%% visualise spread of matched T1s and T2s
% figure; hist(squeeze(matchedT1(offsetListNum,:)))
% figure;
% hist(squeeze(matchedT2(offsetListNum,:)))
% %
% figure; plot(squeeze(data(1,1,1,:)), '-*')
% hold on
% plot(squeeze(bestMatch(1,1,:))*data(1,1,1,1),'--.')

%% Once the similarity function has been run, plot and save the T1, T2 and FA deviation maps
for offsetListNum = 2:8
FA_fig = figure; imagesc(squeeze(matchedFAdevInd(:,:,offsetListNum)))
saveas(FA_fig, [workingdir,'/figures/',phantomName,'matchedFAdevInd_offsetList',num2str(offsetListNum),'_phantomName_',phantomName])
saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum),'_phantomName_',phantomName,'.png'])
matlab2tikz([workingdir,'/Users/jallen/Documents/MATLAB/short_project_2/figures/matchedFAdevInd_offsetList',num2str(offsetListNum),'_phantomName_',phantomName])

matchedT1_fig = figure; imagesc(matchedT1(:,:,offsetListNum))
saveas(matchedT1_fig, [workingdir,'/figures/',phantomName,'matchedT1_offsetList',num2str(offsetListNum),'_phantomName_',phantomName])
saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum),'_phantomName_',phantomName,'.png'])
matlab2tikz([workingdir,'/figures/',phantomName,'matchedT1_offsetList',num2str(offsetListNum),'_phantomName_',phantomName])

matchedT2_fig = figure; imagesc(matchedT2(:,:,offsetListNum))
saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum),'_phantomName_',phantomName])
saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum),'_phantomName_',phantomName,'.png'])
matlab2tikz([workingdir,'/figures/',phantomName,'matchedT2_offsetlist',num2str(offsetListNum),'_phantomName_',phantomName])
end

%%
for offsetListNum = 2:8
plotMap(phantomName,'T1',offsetListNum, workingdir,compartmentCenters)
plotMap(phantomName,'T2',offsetListNum, workingdir,compartmentCenters)
plotMap(phantomName,'FAdevInd',offsetListNum, workingdir,compartmentCenters)
end

%%
%choose pixels and plot time courses for the fingerprinting images
plotTCs(FPimages,22:25, 37, 1, 2) % breaks if 2D array of points chosen

