%% Test the fingerprinting framework with a synthetic phantom
% Author: Jack Allen <jack.allen@jesu.ox.ac.uk>
%
% A script to generate signal timecourses from a synthetic phantom and
% match these timecourses to a dictionary. The idea is to do a sort of
% sanity check on the dictionary generation and matching process as a
% whole. Ideally the matched parameter maps should be the same as the
% synthetic phantom.
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
%%

fingerprintLists = read_FP_offset_list(workingdir,27);
save([savingdir,'/MAT-files/fingerprintLists.mat'], 'fingerprintLists')
load('fingerprintLists.mat')

%%
T1 = 282.3;
T2 = 214.8;

freqOffset = 0;
nSlices = 1;
phantomName = 'sphereD170';
nTimeCoursePts = 48;
%offsetList = generate_offset_list([900, 500, 90, 170],[1100 600 100 190],24,1,'SE_EPI');
offsetListnum = 3;
offsetList = squeeze(fingerprintLists(:,:,offsetListnum));

%Make the synthetic phantom
%[SNR] = calc_SNR(squeeze(FPdata(:,:,1));
SNR = 66; % from sphereD170 fp images
[testFPTCs, testT1map, testT2map] = generate_test_FP_data(20,20,nTimeCoursePts,offsetList,freqOffset,'sphere',SNR);

% Test simulation function seperately
%[~, Mxy] = sim_SE(T1, T2, offsetList,(nTimeCoursePts/24),freqOffset, nSlices, 'noPlot');
%testFPTCs(1,1,:) = Mxy;


%Choose the dictionary parameters
paramList = 8;
dictionaryParams = set_dictionary_params(phantomName,paramList);

%Generate the dictionary, using the dictionary parameters you have just
%chosen
[sd,sdelT] = compile_SE_dictionary(offsetList, nTimeCoursePts, nSlices,freqOffset, dictionaryParams);

% Save the dictionary
%save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'dictionary.mat'],'sd')
%save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'signalDictionaryTime.mat'],'sdelT')

%Calculate the simularity between each the 'data' from the synthetic
%phantom and the dictionary entries
% key:[matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime, maxSimilarity, similarity, degeneracyMap] = calc_similarity(data, dictionary, dictionaryParams)
[matchedT1,matchedT2,matchedFAdev, M0fit_grad, bestMatch, matchTime, maxSimilarity, similarity, degeneracyMap] = calc_similarity(testFPTCs,sd,dictionaryParams);

%% Plot the results
plot_TCs(testFPTCs,'fp',bestMatch)

%%
figure
subplot 221
imagesc(testT1map)
caxis([0 300])
c = colorbar('FontSize',15);
c.Label.String = 'ms';
caxis([0 300])
title('Test T1 map','FontSize',18)
axis square

subplot 222
imagesc(testT2map)
c = colorbar('FontSize',15);
c.Label.String = 'ms';
caxis([0 300])
title('Test T2 map','FontSize',18)
axis square

subplot 223
title('Matched T1 map')
imagesc(matchedT1)
c = colorbar('FontSize',15);
c.Label.String = 'ms';
caxis([0 300])
title('Matched T1 map','FontSize',18)
axis square

subplot 224
title('Matched T2 map')
imagesc(matchedT2)
c = colorbar('FontSize',15);
c.Label.String = 'ms';
caxis([0 300])
title('Matched T2 map','FontSize',18)
axis square
