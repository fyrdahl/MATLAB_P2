% Magnetic Resonance Fingerprinting Pipeline
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015

% import fingerprinting data to workspace
%
clear all
%% Load fingerprinting data
% 1. Set up paths
%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/DPhil';
workingdir = '~/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze

savingdir = '~/Documents/DPhil';
addpath(genpath(savingdir));

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

nPts = 96;
date = '20160727';
ID = 'MR_PHYSICS_494';
set = '1';
series = '005';
%protocol = ['pj_ep2d_se_fp_1SLICE_LIST',num2str(offsetListNum)];
%filepath = [savingdir,'/raw_data/NIfTI/',date,'_',ID,'/',set,'_',series,'_',protocol,'_',date,'/',date,'_',series,'_',protocol];
%filepath = [savingdir,'/data/NIfTI/',date,'_',ID,'/1_003_pj_ep2d_se_fp_1SLICE_LIST19_20160727/20160727_003_pj_ep2d_se_fp_1SLICE_LIST19.nii'];
%filepath = [savingdir,'/data/NIfTI/',date,'_',ID,'/1_038_pj_ep2d_se_fp_1SLICE_LIST20_20160727/20160727_038_pj_ep2d_se_fp_1SLICE_LIST20.nii'];
filepath = [savingdir,'/data/NIfTI/',date,'_',ID,'/1_002_pj_ep2d_se_fp_1SLICE_LIST18_20160727/20160727_002_pj_ep2d_se_fp_1SLICE_LIST18.nii'];
FPimages = read_avw(filepath);
nSlices = size(FPimages,3);
sliceNum = 1;
FPdata = squeeze(FPimages(:,:,sliceNum,1:nPts));

%% 3. Save the list of timing offsets used for acquisition
%fingerprintLists = read_FP_offset_list(workingdir,27);
%offsetList = fingerprintLists(:,:,offsetListNum);
offsets = generate_offset_list([0 0 75 165],[50 50 105 195], nPts, [1 2 3 4],'SSEPI','noSave');
%% Test the simulator
% sphereD170 properties measured by preceeding fits
% T1 = 282.3;
% T2 = 214.8;
% df = 0;
% nRuns = 1; % How many cycles through the offset list?
% row = repmat([1:64],64,1)';
% col = repmat([1:64]',1,64)';
% %coords = [row(:)';col(:)'];
% coords = [ 20:30,20:30];
% plot_sim_comparison(FPdata, T1, T2, offsets,1,df)

%% Create dictionary
% specify the list of dictionary parameters
% dictionaryParamList = 6;
% dictionaryParams = set_dictionary_params('sphereD170',dictionaryParamList);
clear dictionaryParams
dictionaryParams(1,:) = [700:5:900,1100:5:1700,2700:5:3300];
dictionaryParams(2,1:67) = [20:5:350];
%dictionaryParams(3,1:21) = 0.8:0.02:1.2; 
dictionaryParams(3,1) = 1;

df = 0; % frequency offset
nRuns = 1; % How many cycles through the offset list?
nPts = nRuns*size(offsets,1); % How many timepoints (e.g. 48 timecourse points for 2 cycles through a list of 24 entries)
[dictionaryOL18PLSBBern,lut] = compile_SE_dictionary_Bernstein(offsets, 130, 32, nRuns, df, dictionaryParams);
%%

FPdata = reshape(FPdata,size(FPdata,1)*size(FPdata,2),size(FPdata,3))';
matchout = templatematch(dictionaryOL18PLSBBern,FPdata);

% Check similarity scores and use dictionary to assign values of T1, T2 etc.
%[matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime, maxSimilarity, degeneracyMap] = calc_similarity(FPdata, dictionaryOL19PLSBBern, dictionaryParams);

%%

T1matched = reshape(lut(1,matchout),64,64);
T2matched = reshape(lut(2,matchout),64,64);
figure, imagesc(T1matched);
figure, imagesc(T2matched);

bestmatches = dictionaryOL18PLSBBern(matchout,:)

for n = 1:4096
   M0(n)= squeeze(FPdata(:,n)')/squeeze(bestmatches(n,:));
end
M0 = reshape(M0,64,64);
figure, imagesc(M0)

%%
%plot_Matched_maps(matchedT1OL19PLSBBern, matchedT2OL19PLSBBern, matchedFAdevOL18PL9, M0fit_gradOL18PL9Bern)


%%
%plot_TCs(FPdata,'fp',bestMatchOL18PL9Bern)
