%% Magnetic Resonance Fingerprinting (MRF) Pipeline
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015


%% Load fingerprinting data
% 1. Set up paths
workingdir = '~/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze
savingdir = '~/Documents/DPhil';
addpath(genpath(savingdir));
plotFlag = 'plot';

% If working on jalapeno00, uncomment the following lines:
% addpath(genpath('/Applications/fsl/'))
% addpath(genpath('/usr/local/fsl/bin'))
% addpath(genpath('/opt/fmrib/fsl/etc/matlab'))

%% 2. Load images
% Choose the offset list to use. List 2 is the original 'random' list of
% offsets. Lists 3:8 are variations on List2, as described below.
%
% 3: 1st TR offset = 10000
% 4: flip Angles = 90 and 180
% 5: TE offset = 20
% 6: TR offset = 1500, TE offset = 20, FA1 = 90
% 7: TR offset = 1500, TE offset = 20
% 8: TR offset = 15000
%
date = '20160426';
offsetListNum = 18;
ID = 'MR_PHYSICS_459';
set = '1';
series = '019';
protocol = ['pj_ep2d_se_fp_1SLICE_LIST',num2str(offsetListNum)];
filepath = [savingdir,'/data/NIfTI/',date,'_',ID,'/',set,'_',series,'_',protocol,'_',date,'/',date,'_',series,'_',protocol];
%
FPimages = read_avw(filepath);
nSlices = size(FPimages,3);
sliceNum = 1;
FPdata = squeeze(FPimages(:,:,sliceNum,:));

%% 3. Save the list of timing offsets used for acquisition
fingerprintLists = read_FP_offset_list(workingdir,22);
offsetList = fingerprintLists(:,:,12);

if strcmp(plotFlag,'noPlot') ==0
    plot_FP_parameters(offsetList);
end
%% Test the signal simulator
% sphereD170 properties measured by 'gold standard' fits
T1 = 282.3;
T2 = 214.8;
df = 0;
row = repmat(1:64,64,1)';
col = repmat(1:64,64,1);
coords = [ 20:30,20:30];
% user drags a rectangle over the pop-up image to select an area of pixels
% to compare.
nRuns = 2; % How many cycles through the offset list?

if strcmp(plotFlag,'noPlot') ==0
    %function plot_sim_comparison(data, T1, T2, offsets,df)
    plot_sim_comparison_Hargreaves(FPdata, T1, T2, offsetList,nRuns,df);
end
%% Specify the list of dictionary parameters
dictionaryParamListNum = 8;
dictionaryParams = set_dictionary_params('sphereD170',dictionaryParamListNum);
%% Create dictionary
df = 0;
nTimeCoursePts = nRuns*size(offsetList,1); % How many timepoints (e.g. 48 timecourse points for 2 cycles through a list of 24 entries)
[dictionaryOL18PL6BernMRPHYSICS459,dictionaryOL18PL6BernMRPHYSICS459_elapsedTime] = compile_SE_dictionary_Hargreaves(offsetList, 130, 32, nRuns, df, dictionaryParams);

%% Matching
% Check similarity scores and use dictionary to assign values of T1, T2 etc.
[matchedT1OL18PL6BernMRPHYSICS459, matchedTOL18PL6BernMRPHYSICS459, matchedFAdevOL18PL6BernMRPHYSICS459, M0fit_gradOL18PL6BernMRPHYSICS459, bestMatchOL18PL6BernMRPHYSICS459, matchTimeOL18PL6BernMRPHYSICS459, maxSimilarity, similarity, degeneracyMap] = ...
    calc_similarity(FPdata, dictionaryOL18PL6BernMRPHYSICS459, dictionaryParams);

%% Plot Matching Results
limits(1,1:2) = [min(dictionaryParams(1,:)), max(dictionaryParams(1,:))];
limits(2,1:2) = [min(dictionaryParams(2,:)), max(dictionaryParams(2,:))];
limits(3,1:2) = [min(dictionaryParams(3,:)), max(dictionaryParams(3,:))];
if strcmp(plotFlag,'noPlot') ==0
    plot_Matched_maps(matchedT1OL18PL6BernMRPHYSICS459, matchedTOL18PL6BernMRPHYSICS459, matchedFAdevOL18PL6BernMRPHYSICS459, M0fit_gradOL18PL6BernMRPHYSICS459,limits)
end
%% Plot FP data against best matches
if strcmp(plotFlag,'noPlot') ==0
    plot_TCs(FPdata,'fp',bestMatchOL18PL9Bern)
end