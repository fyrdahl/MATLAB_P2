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

date = '2016072';
offsetListNum = 18;
ID = 'MR_PHYSICS_444';
set = '1';
series = '005';
protocol = ['pj_ep2d_se_fp_1SLICE_LIST',num2str(offsetListNum)];
filepath = [savingdir,'/raw_data/NIfTI/',date,'_',ID,'/',set,'_',series,'_',protocol,'_',date,'/',date,'_',series,'_',protocol];
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
df = 0;
nRepeats = 2; % How many cycles through the offset list?
row = repmat([1:64],64,1)';
col = repmat([1:64]',1,64)';
%coords = [row(:)';col(:)'];
coords = [ 20:30,20:30];
plot_sim_comparison(FPdata,coords, T1, T2, offsetList,nSlices,nRepeats,df)

%% Create dictionary
% specify the list of dictionary parameters
dictionaryParamList = 6;
dictionaryParams = set_dictionary_params('sphereD170',dictionaryParamList);

%%
df = 0;
nRepeats = 2; % How many cycles through the offset list?
nTimeCoursePts = nRepeats*size(offsetList,1); % How many timepoints (e.g. 48 timecourse points for 2 cycles through a list of 24 entries)
[dictionaryOL18PL9Bern,sdelT] = compile_SE_dictionary_Bernstein(offsetList, nTimeCoursePts, nSlices, df, dictionaryParams);

% Check similarity scores and use dictionary to assign values of T1, T2 etc.
[matchedT1OL18PL9Bern, matchedT2OL18PL9Bern, matchedFAdevOL18PL9Bern, M0fit_gradOL18PL9Bern, bestMatchOL18PL9Bern, matchTimeOL18PL9Bern, maxSimilarity, similarity, degeneracyMap] = calc_similarity(FPdata, dictionaryOL18PL9Bern, dictionaryParams);


%%


plot_Matched_maps(matchedT1OL18PL9Bern, matchedT2OL18PL9Bern, matchedFAdevOL18PL9, M0fit_gradOL18PL9Bern)


%%
plot_TCs(FPdata,'fp',bestMatchOL18PL9Bern)
