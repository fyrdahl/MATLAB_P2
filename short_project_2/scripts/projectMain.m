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
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze

% If working on jalapeno00, uncomment the following lines:
% addpath(genpath('/Applications/fsl/'))
% addpath(genpath('/usr/local/fsl/bin'))
% addpath(genpath('/opt/fmrib/fsl/etc/matlab'))

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
sliceNumber = 2;
run('plotSamplePixels_TE_TR.m')
%%
run('visualiseImages.m')

%% fit curves to calculate T1 and T2

[compartmentT1s, compartmentT2s, T2curves, T1curves, fittedCurve, goodness, output, F] = fitEvolutionCurves(phantomName,TEimages, TIimages, TE(2:end)', TI(2:end), 'compartments', compartmentCenters);

%% 6. read in the list of timing offsets used for acquisition
run('readFingerprintOffsetList.m')
save([workingdir,'/MAT-files/fingerprintLists.mat'], 'fingerprintLists')
%% 7. simulate magnetisation evolution
%check bloch simulation by using properties of the phantom
load('fingerprintLists.mat')
% sphereD170 properties
T1 = 282.3;
T2 = 214.8;
% fingerprintLists(:,1,offsetListNum) = 100;
% fingerprintLists(:,2,offsetListNum) = 200;
% fingerprintLists(:,3,offsetListNum) = 90;
% fingerprintLists(:,4,offsetListNum) = 180;
freqOffset = 0;
nSlices = 2;
offsetListNum = 3;
nTimeCoursePts = 24;
[M, Mxy,flipAngles, t0s] = SimBloch(T1, T2, fingerprintLists(:,:,offsetListNum), 'showPlot', freqOffset, nSlices, nTimeCoursePts);
%[M, Mxy,flipAngles, t0s] = SimBloch2(T1, T2, fingerprintLists(:,:,offsetListNum), 'showPlot', freqOffset, nSlices);

% 8. check signal simulation by plotting positions of sample pixels for the fingerprinting images
compartmentCentersList = 3;
clear FPimages
load ([phantomName,'_list',num2str(offsetListNum),'FPimages.mat'])
run('plotSamplePixels.m')

% plot comparison of simulation with sampled pixels
run('plotSim.m')

%% 9. create dictionary
load fingerprintLists.mat
%
paramList = 3
[dictionaryParams, paramList] = setDictionaryParams(phantomName,paramList);
%
nTimeCoursePts = 24;
%%
for offsetListNum = 2:8;
    %%
    [signalDictionary, sdelT] = compileDictionary(fingerprintLists, offsetListNum, dictionaryParams, nTimeCoursePts, freqOffset, nSlices, phantomName, background);
    save([workingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'dictionary.mat'],'signalDictionary')
    save([workingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'signalDictionaryTime.mat'],'sdelT')
    
    pause(1)
end

%% 10. check similarity and use dictionary to measure T1 and T2

sliceNumber = 2 ;% slice to be analysed
clear data

% data = FPimages(compartmentCenters(1,1),compartmentCenters(1,2),:,:);

% for r = 1:size(compartmentCenters,1)
% data(r,1,sliceNumber,:) = FPimages(compartmentCenters(r,1),compartmentCenters(r,2),sliceNumber,:);
% end
%% SIMILARITY
for offsetListNum = 2:8;
    %%
    sliceNumber
    clear FPimages
    clear Mzeq
    clear data
    clear M0_stdPC
    clear M0
    clear M0_mean
    clear matchedT1
    clear similarity
    clear matchedT2
    clear matchedFAdevInd
    load([workingdir,'/MAT-files/images/',phantomName,'_list',num2str(offsetListNum),'FPimages.mat'])
    data = zeros(24, size(FPimages,1)*size(FPimages,2));
    %   data(r,c,sliceNumber,1:24) = squeeze(FPimages(r,c,sliceNumber,1:24,offsetListNum));
    data = squeeze(FPimages(:,:,sliceNumber,1:24,offsetListNum));
    
    disp(['offsetList',num2str(offsetListNum),', phantom:',phantomName])
    clear signalDictionary
    l = load([workingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'dictionary.mat']);
    signalDictionary = l.signalDictionary;
    [similarity, matchedT1, matchedT2, matchedFAdevInd, M0_mean, M0_stdPC,M0,  scales, bestMatch, matchT] = calcSimilarity(data, signalDictionary(:,:,:,:,offsetListNum), sliceNumber, dictionaryParams, workingdir);
   
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedT1.mat'],'matchedT1')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedT2.mat'],'matchedT2')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedFAdevInd.mat'],'matchedFAdevInd')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0_mean.mat'],'M0_mean')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0_stdPC.mat'],'M0_stdPC')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0.mat'],'M0')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'scales.mat'],'scales')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'bestMatch.mat'],'bestMatch')
    save([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'compileDictionaryElapsedTime.mat'],'matchT')
    
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
    %%
    offsetListNum
    clear FA_fig
    clear matchedT1_fig
    clear matchedT2_fig
    clear matchedT1
    clear matchedT2
    clear matchedFAdevInd
    clear M0_mean
    load([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedT1.mat'])
    load([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedT2.mat'])
    load([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedFAdevInd.mat'])
    load([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0_mean.mat'])
    load([workingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0.mat'])
    
    filename = '/Users/jallen/Documents/MATLAB/short_project_2/DTC_report/';
    
    
    FA_fig = figure
    set(FA_fig,'name',[phantomName,', List',num2str(offsetListNum),', B1 deviation'])
    imagesc(squeeze(matchedFAdevInd(:,:)))
    axis off
    colormap jet
    %   title ([phantomName,', List',num2str(offsetListNum),', B1 deviation'])
    colorbar
    cFA = colorbar;
    cmin = min(dictionaryParams(3,1:sum(dictionaryParams(3,:) ~= 0))) - 0.1*min(dictionaryParams(3,1:(sum(dictionaryParams(3,:) ~= 0))));
    cmax = max(dictionaryParams(3,1:sum(dictionaryParams(3,:) ~= 0)))
    %  cFA.YTick = [cmin , cmax]
    caxis([cmin,cmax])
    ylabel(cFA,'Fraction of B1')
    
    %saveas(FA_fig, [workingdir,'/figures/',phantomName,'matchedFAdevInd_offsetList',num2str(offsetListNum),'_phantomName_',phantomName])
    saveas(FA_fig, [filename,'/',phantomName,'matchedFAdevIndoffsetlist',num2str(offsetListNum),'.png'])
    matlab2tikz('figurehandle',FA_fig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'FAlist',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    
    pause(2)
    matchedT1_fig = figure
    set(matchedT1_fig,'name',[phantomName,', List',num2str(offsetListNum),', T1'])
    imagesc(matchedT1(:,:))
    axis off
    colormap jet
    %title ([phantomName,', List',num2str(offsetListNum),', T1'])
    cT1 = colorbar;
    switch phantomName
        case 'sphereD170'
            cmin = min(dictionaryParams(1,1:sum(dictionaryParams(1,:) ~= 0))) - 0.1*min(dictionaryParams(1,1:(sum(dictionaryParams(1,:) ~= 0))));
            cmax = max(dictionaryParams(1,1:sum(dictionaryParams(1,:) ~= 0)))
            
        case 'Jack'
            cmax = compartmentT1s(2)
            cmin = 50
            cmax = 260
    end
    cT1.YTick = [cmin : 10 : cmax]
    caxis([cmin,cmax])
    ylabel(cT1,'T1 (ms)')
    
    %saveas(matchedT1_fig, filenameT1 )
    %saveas(matchedT1_fig, [filename,'/',phantomName,'matchedT1offsetlist',num2str(offsetListNum),'-1.png'])
    matlab2tikz('figurehandle',matchedT1_fig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'T1list',num2str(offsetListNum)','ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    
    clear temp
    pause(2)
    matchedT2_fig = figure; imagesc(matchedT2(:,:))
    axis off
    set(matchedT2_fig,'name',[phantomName,', List',num2str(offsetListNum),', T2'])
    colormap jet
    %title ([phantomName,', List',num2str(offsetListNum),', T2'])
    cT2 = colorbar;
    switch phantomName
        case 'sphereD170'
            cmin = min(dictionaryParams(2,1:sum(dictionaryParams(2,:) ~= 0))) - 0.1*min(dictionaryParams(2,1:(sum(dictionaryParams(2,:) ~= 0))));
            cmax = max(dictionaryParams(2,1:sum(dictionaryParams(2,:) ~= 0))) ;
        case 'Jack'
            for i = 1:size(compartmentCenters(:,:,3),1)-1
                temp(i) = matchedT2(squeeze(compartmentCenters(i,1,2)),squeeze(compartmentCenters(i,2,2)))
                % tp(i) = squeeze(compartmentCenters(i,:,3));
            end
            cmin = min(temp)
            cmax = max(temp)
            cmin = 10
            cmax = 110
            
    end
    %   cmax = 120
    % cmax = compartmentT2s(2) + 0.1*compartmentT2s(2)
    cT2.YTick = [cmin : 10 : cmax];
    caxis([cmin,cmax])
    %  caxis([min(dictionaryParams(2,1:sum(dictionaryParams(2,:) ~= 0))) - 0.1*min(dictionaryParams(2,1:(sum(dictionaryParams(2,:) ~= 0)))) ,max(dictionaryParams(2,1:(sum(dictionaryParams(2,:) ~= 0)))) ])
    ylabel(cT2,'T2 (ms)')
    %saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum)])
    saveas(matchedT2_fig, [filename,'/',phantomName,'matchedT2offsetlist',num2str(offsetListNum),'.png'])
    matlab2tikz('figurehandle',matchedT2_fig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'T2list',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    pause(2)
    
    image = log10(abs(M0_mean(:,:)));
   M0_mean_fig = figure; imagesc(image(:,:))
   axis off
    set(M0_mean_fig,'name',[phantomName,', List',num2str(offsetListNum),', M0_mean'])
    colormap jet
    cM0_mean = colorbar;
             %    for i = 1:size(compartmentCenters(:,:,3),1)
             %        temp(i) = image(squeeze(compartmentCenters(i,1,3)),squeeze(compartmentCenters(i,2,3)));
             %         tp(i) = squeeze(compartmentCenters(i,:,3));
             %   end
             %   cmin = min(temp)
             %   cmax = max(temp)
                %cmax = 120
             %   caxis([cmin,cmax])
   ylabel(cM0_mean,'M0 (log_{10}(Mean of Absolute Scaling Factors))')
    saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum)])
   saveas(M0_mean_fig, [filename,'/',phantomName,'matchedM0_mean',num2str(offsetListNum),'.png'])
   matlab2tikz('figurehandle',M0_mean_fig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'M0mean',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    
    %     pause(2)
    %     M0fig = figure; imagesc(log10(M0(:,:)))
    %     axis off
    %     set(M0fig,'name',[phantomName,', List',num2str(offsetListNum),', M0'])
    %     colormap jet
    %     cT2 = colorbar;
    %     ylabel(cT2,'log_{10}(M0)')
    %
    %     %saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum)])
    %     saveas(M0fig, [filename,'/',phantomName,'matchedM0',num2str(offsetListNum),'.png'])
    %     matlab2tikz('figurehandle',M0fig,'filename',[filename,'/',phantomName,'M0',num2str(offsetListNum)],'height', '\figureheight', 'width', '\figurewidth')
    %
end

%%
for offsetListNum = 2:8
    
    plotMap(phantomName,'T1',offsetListNum, workingdir,compartmentCenters)
    plotMap(phantomName,'T2',offsetListNum, workingdir,compartmentCenters)
    plotMap(phantomName,'FAdevInd',offsetListNum, workingdir,compartmentCenters)
end

%%
%choose pixels and plot time courses for the fingerprinting images
plotTCs(FPimages,15:4:30, 37, 1, 2) % breaks if 2D array of points chosen

%%
for offsetListNum = 2:8
    load(['/Users/jallen/Documents/MATLAB/short_project_2/MAT-files/matches/Jacklist',num2str(offsetListNum),'paramList1scales.mat'])
    scalesFig = figure;
    for i = 1:size(compartmentCenters(:,1))
        plot(log10(squeeze(scales(compartmentCenters(i,1,3),compartmentCenters(i,2,3),:))),'.-')
        hold on
    end
    ylabel (['log_{10}(Scaling Factor)'])
    xlabel (['Image Index'])
    legend ({'Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Compartment 5', 'Compartment 6'},'Position',[0.35,0.6,0.25,0.1],'FontSize',8)
    filename = '/Users/jallen/Documents/MATLAB/short_project_2/DTC_report/';
    
    matlab2tikz('figurehandle',scalesFig,'filename',[filename,'/',phantomName,'compartmentScales',num2str(offsetListNum)],'height', '\figureheight', 'width', '\figurewidth')
    
end