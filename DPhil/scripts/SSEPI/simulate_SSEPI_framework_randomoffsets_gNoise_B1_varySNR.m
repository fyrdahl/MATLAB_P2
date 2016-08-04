%% Script to simulate the performance of SSEPI fingerprinting sequence
%
% Jack Allen <jack.allen@jesus.ox.ac.uk>
clear all
disp('Simulate Performance of SSEPI MRF sequence...')
%% make synthetic data
T1 = 282;
T2 = 214;
dim = 30;
T1 = phantom('Modified Shepp-Logan',dim)*T1;
T2 = phantom('Modified Shepp-Logan',dim)*T2;
%
% B1 distribution
B1 = fspecial('gaussian', dim, dim/2);
B1 = B1./(max(max(B1)));

%% generate offset list
nPts = 10;
% add hard-coded sequence timings
minTR = 130; % ms, from sequence protocol card.
minTE = 32; % ms, from sequence protocol card.
%   offsets = generate_offset_list(lower,upper,nPts,seed,seq,saveFlag);
offsets = generate_offset_list([0 0 75 165],[50 50 105 195],nPts,[1 2 3 4],'SSEPI','noSave');

%% make dictionary
% Specify the list of dictionary parameters
T1s = 0:5:300;
T2s = 0:5:300;
FAdevs = 1;% 0.75:0.01:1.25;
dictionaryParams(1,1:numel(T1s)) = T1s; % T1
dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction

nRuns = 1;
df = 0;
[dictionary, paramIndices] = compile_SE_dictionary_Bernstein(offsets, minTR, minTE, nRuns, df, dictionaryParams);

% simulate acquired data
% [testFPtcMap, testT1map, testT2map] = generate_test_FP_data(dim1, dim2,nTimePts,offsetList,freqOffset,mapType,SNR)
SNRs = [50:50:200];
for SNR = SNRs;
    
    nSNR = find(SNR == SNR);
syntheticFPdata = generate_test_FP_data(T1,T2, minTR, minTE, nRuns,offsets,df,SNR,B1);

matchout=templatematch(dictionary,syntheticFPdata);
matchedT1 = paramIndices(1,matchout);
matchedT2 = paramIndices(2,matchout);
matchedB1 = paramIndices(3,matchout);


% RMSE 
for nVoxel = 1:numel(matchedT1)
T1err(nVoxel) = matchedT1(nVoxel) - T1(nVoxel);   % Errors
T2err(nVoxel) = matchedT2(nVoxel) - T2(nVoxel);   % Errors
end
T1RMSE = sqrt(mean((T1err).^2));  % Root Mean Squared Error
T2RMSE = sqrt(mean((T2err).^2));  % Root Mean Squared Error

allT1RMSE(nSNR,:) = T1RMSE;
allT2RMSE(nSNR,:) = T2RMSE;
%
end
matchedT1 = reshape(matchedT1,dim,dim);
matchedT2 = reshape(matchedT2,dim,dim);
matchedB1 = reshape(matchedB1,dim,dim);

%% simulate acquired data
% [testFPtcMap, testT1map, testT2map] = generate_test_FP_data(dim1, dim2,nTimePts,offsetList,freqOffset,mapType,SNR)
SNR = 75;
syntheticFPdata = generate_test_FP_data(T1,T2, minTR, minTE, nRuns,offsets,df,SNR,B1);


%% Compare dictionary with simulated signal
addpath '~/Documents/MATLAB/Ma2013'
matchout=templatematch(dictionary,syntheticFPdata);
matchedT1 = paramIndices(1,matchout);
matchedT2 = paramIndices(2,matchout);
matchedB1 = paramIndices(3,matchout);


% RMSE 
for nVoxel = 1:numel(matchedT1)
T1err(nVoxel) = matchedT1(nVoxel) - T1(nVoxel);   % Errors
T2err(nVoxel) = matchedT2(nVoxel) - T2(nVoxel);   % Errors
end
T1RMSE = sqrt(mean((T1err).^2));  % Root Mean Squared Error
T2RMSE = sqrt(mean((T2err).^2));  % Root Mean Squared Error

%%
matchedT1 = reshape(matchedT1,dim,dim);
matchedT2 = reshape(matchedT2,dim,dim);
matchedB1 = reshape(matchedB1,dim,dim);


%% Plot Results
disp 'Plotting Results...'
T1 = reshape(T1,dim,dim);
T2 = reshape(T2,dim,dim);
matchedT1 = reshape(matchedT1,dim,dim);
matchedT2 = reshape(matchedT2,dim,dim);
%
%
figure('name',['SSEPI - with Gaussian Noise (SNR = ',num2str(SNR),' nPts = ',num2str(nPts),')'])
%
subplot 341
imagesc(T1);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Synthetic T1'
%
subplot 345
imagesc(T2);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Synthetic T2'
%
subplot 342
imagesc(matchedT1);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Matched T1'
%
subplot 346
imagesc(matchedT2);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Matched T2'
%

subplot 343
imagesc(abs(matchedT1-T1))
%caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title(['Error (RMSE = ',num2str(T1RMSE),')'])
%
subplot 347
imagesc(abs(matchedT2-T2))
%caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title(['Error (RMSE = ',num2str(T2RMSE),')'])
%
subplot 344
plot(T1(:),matchedT1(:),'o')
hold on
plot(0:300,0:300)
xlim([0 300])
ylim([0 300])
ylabel(c,'T1 [ms]')
ylabel 'Matched T1'
xlabel 'Synthetic T1'
%
subplot 348
plot(T2(:),matchedT2(:),'o')
hold on
plot(0:300,0:300)
xlim([0 300])
ylim([0 300])
ylabel(c,'T2 [ms]')
ylabel 'Matched T2'
xlabel 'Synthetic T2'
%
subplot 349
imagesc(B1);
c=colorbar;
colormap hot
ylabel(c,'Efficiency')
title 'Synthetic B1'

disp 'Plotting Results... Done.'
