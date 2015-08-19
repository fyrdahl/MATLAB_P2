clear all
close all
addpath(genpath('/Applications/fsl/'))
addpath(genpath('/usr/local/fsl'))
addpath(genpath('/usr/local/fsl/bin'))
addpath(genpath('/Users/jallen/Documents/MATLAB'))
%%

%which phantom is the data from? ('sphereD170' or 'Jack'?)
phantomName = 'Jack'
[TEImageInfo, TIImageInfo, FPImageInfo, TEimages, TIimages, FPimages, TE, TI] = readData(phantomName, 'offsetList2');
%%
[SNR signal noise] = calcSNR(TEimages,TE,'showFigure');
%%
run('compartmentSignals.m')
%%
% mask including whole phantom
load('/Users/jallen/Documents/MATLAB/short_project_2/mask.mat')
run('visualiseImages.m')

%% normalise pixel values
% TEimages(:,:,TE) = TEimages(:,:,TE)/max(max(max(TEimages(:,:,:))));
% TIimages(:,:,TE) = TIimages(:,:,TI)/max(max(max(TIimages(:,:,:))));

%% Calculate T1 and T2

[ fittedCurve, goodness, output] = fitEvolutionCurves(TEimages, 'T2', TE(2:end), 'compartments', compartmentCenters)

[ fittedCurve, goodness, output] = fitEvolutionCurves(TIimages, 'T1', TI(2:end), 'compartments', compartmentCenters)

%% Calculate T1 and T2: with means

% set the position and size of the ROI from which to calculate the mean
ROI_dim = [20:25; 20:25];
ROI_dim = [23:26; 27:30];

for i = TE
ROI = (TEimages(ROI_dim(1,:),ROI_dim(2,:),i));
TEmeans(i) = (mean(ROI(:)));
end
for i = TI
ROI = (TIimages(ROI_dim(1,:),ROI_dim(2,:),i));
TImeans(i) = (mean(ROI(:)));
end
run('calcT1T2withMeans.m')

%%
run('readFingerprintOffsetList.m')

%% simulate magnetisation evolution

%check bloch simulation by using properties of the phantom
% 
% fingerprintList1(:,1) = 3000
% fingerprintList1(:,2) = [100]
% fingerprintList1(:,3) = 90
% fingerprintList1(:,4) = 180

nSlices = 2;
T1 = 275.26;
T2 = 214.2;

freqOffset = 0
[simMtransverse, simImageMtransverse, M,TE,TR] = SimBloch(T1, T2, fingerprintList2, 'showPlot',freqOffset);
%% 
simSignal(2,1:(size(fingerprintList2(:,1),1))) = simImageMtransverse;
%%
nSlices = 2;
[M,simSignal,TE,TR,imageTimes,newTRs] = simulateSignal(T1, T2, nSlices, fingerprintList2, 'showPlots');

%%
figure; hold on
 ySim = simImageMtransverse(1:end)./simImageMtransverse(1);
% plot(ySim,'o')
% ySim = simSignal(2,:)/simSignal(2,1);
plot(ySim,'x','MarkerSize',20)
y = zeros(3,(size(FPimages,4)/2));
sliceNumber = 1;

ROIimageMean = mean(mean(squeeze(FPimages(40:45,20:40,sliceNumber,:))));
ROIimageMean(:) = ROIimageMean/ROIimageMean(1);
yImageROI = squeeze(ROIimageMean);
plot(yImageROI(1:24),'^')
for n = 1:6
for i = 1:(size(FPimages,4)/2)
y(n,i) = squeeze(FPimages(compartmentCenters(n,1),compartmentCenters(n,2),sliceNumber,i));
end

y(n,:) = y(n,:)/y(n,1);
% residuals = y(n,:) - ySim;
plot(y(n,1:(size(FPimages,4)/2)),'*')
end
legend 'Simulated Signal'
xlabel 'TE indices'
ylabel 'signal (normalised to first measurement)'
%  figure; plot(residuals,'+')

%% simularity measure
for n = 1:6
% Cosine similarity (does not depend on magnitude of each vector)
similarity(n) = dot(ySim,y(n,:))/(norm(ySim)*norm(y(n,:)))
end


%% DICTIONARY
T1dictionaryList = [50;100;150;200;250;300]
T2dictionaryList = [50;100;150;200;250;300]
for T1 = T1dictionaryList'
    T1
    for T2 = T2dictionaryList'  
        T2
[Mtransverse(find(T1dictionaryList == T1),find(T2dictionaryList == T2)) M] = SimBloch(T1, T2, fingerprintList2, 'showPlot'); 
    end
end





