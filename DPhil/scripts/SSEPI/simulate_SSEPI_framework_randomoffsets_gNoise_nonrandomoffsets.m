%% Script to simulate the performance of SSEPI fingerprinting sequence
%
% Jack Allen <jack.allen@jesus.ox.ac.uk>
clear all
disp('Simulate Performance of SSEPI MRF sequence...')
%% make synthetic data
T1 = 282;
T2 = 214;
dim = 64;
T1 = phantom('Modified Shepp-Logan',dim)*T1;
T2 = phantom('Modified Shepp-Logan',dim)*T2;
T1 = awgn(T1,10);
T2 = awgn(T2,10);
T1(T1<1) =0;
T2(T2<1) =0;

%% generate offset list
nPts = 48;
% add hard-coded sequence timings
minTR = 130; % ms, from sequence protocol card.
minTE = 32; % ms, from sequence protocol card.

workingdir = '/Users/jallen/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze
savingdir = '/Users/jallen/Documents/DPhil';
addpath(genpath(savingdir));
fingerprintLists = read_FP_offset_list(workingdir,27);
offsets = fingerprintLists(:,:,2);

%% make dictionary
% Specify the list of dictionary parameters
T1s = 0:5:300;
T2s = 0:5:300;
FAdevs =  1;
dictionaryParams(1,1:numel(T1s)) = T1s; % T1
dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction

nRuns = nPts/size(offsets,1);
df = 0;
[dictionary, paramIndices] = compile_SE_dictionary_Bernstein(offsets, minTR, minTE, nRuns, df, dictionaryParams);

%% simulate acquired data
% [testFPtcMap, testT1map, testT2map] = generate_test_FP_data(dim1, dim2,nTimePts,offsetList,freqOffset,mapType,SNR)
SNR = 75;
syntheticFPdata = generate_test_FP_data(T1,T2, minTR, minTE, nRuns,offsets,df,SNR);


%% Compare dictionary with simulated signal
addpath '~/Documents/MATLAB/Ma2013'
matchout=templatematch(dictionary,syntheticFPdata);
matchedT1 = paramIndices(1,matchout);
matchedT2 = paramIndices(2,matchout);
matchedB1 = paramIndices(3,matchout);

matchedT1 = reshape(matchedT1,dim,dim);
matchedT2 = reshape(matchedT2,dim,dim);
matchedB1 = reshape(matchedB1,dim,dim);


% RMSE 
for nVoxel = 1:numel(matchedT1)
T1err(nVoxel) = matchedT1(nVoxel) - T1(nVoxel);   % Errors
T2err(nVoxel) = matchedT2(nVoxel) - T2(nVoxel);   % Errors
end
T1RMSE = sqrt(mean((T1err).^2));  % Root Mean Squared Error
T2RMSE = sqrt(mean((T2err).^2));  % Root Mean Squared Error



figure('name',['SSEPI - with Gaussian Noise (SNR = ',num2str(SNR),' nPts = ',num2str(nPts),')'])
%
subplot 241
imagesc(T1);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Synthetic T1'
axis square
%
subplot 245
imagesc(T2);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Synthetic T2'
axis square

%
subplot 242
imagesc(matchedT1);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Matched T1'
axis square

%
subplot 246
imagesc(matchedT2);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Matched T2'
axis square

%

subplot 243
imagesc(abs(matchedT1-T1))
%caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title(['Error (RMSE = ',num2str(T1RMSE),')'])
axis square

%
subplot 247
imagesc(abs(matchedT2-T2))
%caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title(['Error (RMSE = ',num2str(T2RMSE),')']) 
axis square

%
subplot 244
plot(T1(:),matchedT1(:),'o')
hold on
plot(0:300,0:300)
xlim([0 300])
ylim([0 300])
ylabel(c,'T1 [ms]')
ylabel 'Matched T1'
xlabel 'Synthetic T1'
axis square

%
subplot 248
plot(T2(:),matchedT2(:),'o')
hold on
plot(0:300,0:300)
xlim([0 300])
ylim([0 300])
ylabel(c,'T2 [ms]')
ylabel 'Matched T2'
xlabel 'Synthetic T2'
axis square

%

%
disp 'Plotting Results... Done.'


%%
syntheticFPdata = reshape(syntheticFPdata,nPts,dim,dim);
matchout = reshape(matchout,dim,dim);
fig = figure;
subplot 121
title ('example image (image 1)')
imagesc(squeeze(syntheticFPdata(1,:,:)));
[row, col] = getpts(fig);
row = round(row);
col = round(col);

figure
for n= 1:numel(row)
plot(squeeze(syntheticFPdata(:,col(n),row(n))))
hold on
plot(squeeze(dictionary(matchout(col(n),row(n)),:)),'-')
end