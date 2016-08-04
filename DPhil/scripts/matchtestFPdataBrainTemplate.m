%% Test the fingerprinting framework with a synthetic phantom
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
%
% A script to generate signal timecourses from a synthetic phantom and
% match these timecourses to a dictionary. The idea is to do a sort of
% sanity check on the dictionary generation and matching process as a
% whole. Ideally the matched parameter maps should be the same as the
% synthetic phantom.


freqOffset = 0;
nSlices = 1;
phantomName = 'brainTemplate';
nTimeCoursePts = 48;
paramList = 'brainTemplate';

% Generate a list of offsets for each parameter in the sequence (e.g. TE,
% TR etc)
offsetList = generate_offset_list([0, 0, 90, 170],[800 600 100 190],24,1,'SE_EPI');

%Make the synthetic phantom
[testFPTCs, testT1map, testT2map] = generate_test_FP_data(nTimeCoursePts,offsetList,freqOffset,'brainTemplate');

% Test simulation function seperately
%[~, Mxy] = sim_SE(T1, T2, offsetList,(nTimeCoursePts/24),freqOffset, nSlices, 'noPlot');
%testFPTCs(1,1,:) = Mxy;

%Choose the dictionary parameters
[dictionaryParams] = set_dictionary_params(phantomName);

%Generate the dictionary, using the dictionary parameters you have just
%chosen
[sd,sdelT] = compile_SE_dictionary(offsetList, nTimeCoursePts, freqOffset, dictionaryParams);

% Save the dictionary
%save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'dictionary.mat'],'sd')
%save([savingdir,'/MAT-files/dictionaries/',phantomName,'_list',num2str(offsetListNum),'paramList',num2str(paramList),'signalDictionaryTime.mat'],'sdelT')

%Calculate the simularity between each the 'data' from the synthetic
%phantom and the dictionary entries
% [maxSimilarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, matchTime] = calc_similarity(data, sd,nTimeCoursePts, paramRangeList,savingdir, inputPhantomName,offsetListNum)
[similarity, matchedT1, matchedT2, matchedFAdev, M0fit_grad, bestMatch, match_time] = calc_similarity(testFPTCs,sd,dictionaryParams);

%% Plot the results
plot_TCs(testFPTCs,'fp',bestMatch)

%%
            
figure
subplot 221
imagesc(testT1map)
caxis([700 3200])
colorbar
title('Test T1 map')

subplot 222
imagesc(testT2map)
caxis([0 600])
colorbar
title('Test T2 map')

subplot 223
title('Matched T1 map')
imagesc(matchedT1)
colorbar
caxis([700 3200])
title('Matched T1 map')

subplot 224
title('Matched T2 map')
imagesc(matchedT2)
caxis([0 600])
colorbar
title('Matched T2 map')