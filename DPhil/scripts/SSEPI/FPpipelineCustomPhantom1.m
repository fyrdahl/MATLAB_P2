%% Magnetic Resonance Fingerprinting (MRF) Pipeline for Custom Phantom 1 ('Jack1')
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015

%% 1. Set up paths
%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/DPhil';
workingdir = '/Users/jallen/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze

savingdir = '/Users/jallen/Documents/DPhil';
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

ID = 'MR_PHYSICS_352';
set = '2';
series = '072';
offsetListNum = 2;
protocol = ['pj_ep2d_se_fp_List',num2str(offsetListNum),'_2sl_tr132start'];
date = '20150806';
%which slice should we analyse?
sliceNum = 1; 

FPimages = nifti2mat([savingdir,'/raw_data/NIfTI/',date,'_',ID,'/',set,'_',series,'_',protocol,'_',date,'/',date,'_',series,'_',protocol]);
FPdata = squeeze(FPimages(:,:,sliceNum,:));

%% 3. Save the list of timing offsets used for acquisition
fingerprintLists = read_FP_offset_list(workingdir);
save([savingdir,'/MAT-files/fingerprintLists.mat'], 'fingerprintLists')
load('fingerprintLists.mat')
offsetList = fingerprintLists(:,:,offsetListNum);
%% 4. Settings
nTimeCoursePts = nRepeats*size(fingerprintLists,1); % How many timepoints
freqOffset = 0;

coords = [23; 23];
plot_sim_comparison(data,coords, T1, T2,offsetList,nRepeats,0)

%% Create dictionary
% specify the list of dictionary parameters
dictionaryParamListNum = 4;
dictionaryParams = set_dictionary_params('Jack1',dictionaryParamListNum);
[dictionary,sdelT] = compile_SE_dictionary(offsetList, nTimeCoursePts, freqOffset, dictionaryParams);

%% Check similarity scores and use dictionary to assign values of T1, T2 etc.
% [maxSimilarity, degeneracyMap, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime] = calc_similarity(data, sd, dictionaryParams)
[maxSimilarity,degeneracyMap, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, match_time] = calc_similarity(data,dictionary, dictionaryParams);
