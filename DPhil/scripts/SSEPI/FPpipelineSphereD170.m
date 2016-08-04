% Magnetic Resonance Fingerprinting Pipeline
% Author: Jack Allen
% Supervisor: Prof. Peter Jezzard
% Start Date: 13th July 2015

% import fingerprinting data to workspace
%

%% Load fingerprinting data
% 1. Set up paths
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

date = '20150819';
offsetListNum = 4;
ID = 'MR_PHYSICS_358';
set = '1';
series = '00';
protocol = ['pj_ep2d_se_fp_List',num2str(offsetListNum)];
filepath = [savingdir,'/data/NIfTI/',date,'_',ID,'/',set,'_',series,num2str(offsetListNum+1),'_',protocol,'_',date,'/',date,'_',series,num2str(offsetListNum+1),'_',protocol,'.nii'];
FPimages = read_avw(filepath);
nSlices = size(FPimages,3);
sliceNum = 1;
FPdata = squeeze(FPimages(:,:,sliceNum,:));

%% 3. Save the list of timing offsets used for acquisition
fingerprintLists = read_FP_offset_list(workingdir,27);
offsetList = fingerprintLists(:,:,offsetListNum);

%% Test the simulator
% sphereD170 properties measured by preceeding fits

T1 = 282.3;
T2 = 214.8;
df = 50;
nRepeats = 2; % How many cycles through the offset list?
plot_sim_comparison(FPdata, T1, T2, offsetList,nSlices,df)

%% Create dictionary
% specify the list of dictionary parameters
dictionaryParamList = 8;
dictionaryParams = set_dictionary_params('sphereD170',dictionaryParamList);

%%
df = 0;
nRepeats = 2; % How many cycles through the offset list?
nTimeCoursePts = nRepeats*size(offsetList,1); % How many timepoints
[dictionaryOL3PL8,sdelT] = compile_SE_dictionary(offsetList, nTimeCoursePts, nSlices, df, dictionaryParams);

% Check similarity scores and use dictionary to assign values of T1, T2 etc.
[maxSimilarity, similarity,degeneracyMap, matchedT1OL3PL8, matchedT2OL3PL8, matchedFAdevOL3PL8, M0fit_gradOL3PL8, bestMatch, matchTime] = calc_similarity(FPdata, dictionary, dictionaryParams);

%%
plot_TCs(FPdata,'fp',bestMatch)
