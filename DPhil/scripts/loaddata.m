%% loaddata
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
offsetListNum = 3;

%which slice should we analyse?
sliceNumber = 2; 

readData(phantomName, offsetListNum, [savingdir,'/raw_data/'], savingdir );
[TEimages, TIimages, TE, TI] = load_GS_Data(phantomName, offsetListNum, savingdir);
[FPimages] = load_FP_Data(phantomName, offsetListNum, sliceNumber, savingdir);

%% 3. Save the list of timing offsets used for acquisition
fingerprintLists = read_FP_offset_list(workingdir);
save([savingdir,'/MAT-files/fingerprintLists.mat'], 'fingerprintLists')