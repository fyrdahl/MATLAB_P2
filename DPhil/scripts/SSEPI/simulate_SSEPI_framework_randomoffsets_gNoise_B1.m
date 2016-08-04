%% Script to simulate the performance of SSEPI fingerprinting sequence
%
% Jack Allen <jack.allen@jesus.ox.ac.uk>
clear all
disp('Simulate Performance of SSEPI MRF sequence...')
%% Make Synthetic Phantom
SNR = 75;
addpath('~/Documents/MATLAB/Ma2013')
dim = 32;
refT1 = 282;
refT2 = 214;

phantom = 'SL';

switch phantom
    
    case 'SL'
        if dim > 1;
            T1 = round(phantom('Modified Shepp-Logan',dim)*10);
            T2 = round(phantom('Modified Shepp-Logan',dim)*10);
            T1(T1<2) = 0;
            T2(T2<2) = 0;
            
            T1(T1==2) = 1*refT1;
            T1(T1==3) = 0.9*refT1;
            T1(T1==4) = 0.95*refT1;
            T1(T1==10) = 0.9*refT1;
            
            T2(T2==2) = 1*refT2;
            T2(T2==3) = 1.05*refT2;
            T2(T2==4) = 1*refT2;
            T2(T2==10) = 1.1*refT2;
            
            T1 = awgn(T1,5);
            T2 = awgn(T2,5);
            
            T1(T1<10) = 0;
            T2(T2<10) = 0;
        end
        
    case 'Brain' %Brain Values
        T1 = ones(dim);
        T2 = ones(dim);
        
        [wmT1, wmT2] = get_relaxation_times(3,'wm');
        [gmT1, gmT2] = get_relaxation_times(3,'gm');
        [bloodT1, bloodT2] = get_relaxation_times(3,'blood');
        csfT1 = 3000;
        csfT2 = 300;
        T1(:,1:dim/4) = wmT1;
        T1(:,dim/4:dim/2) = gmT1;
        T1(:,dim/2:(3*dim/4)) = bloodT1;
        T1(:,(3*dim/4):end) = csfT1;
        
        T2(:,1:dim/4) = wmT2;
        T2(:,dim/4:dim/2) = gmT2;
        T2(:,dim/2:(3*dim/4)) = bloodT2;
        T2(:,(3*dim/4):end) = csfT2;
        
    case 'simpleBrain' %Brain Values
        
        nTissue = 4;
        T1 = ones(1,nTissue);
        T2 = ones(1,nTissue);
        
        [wmT1, wmT2] = get_relaxation_times(3,'wm');
        [gmT1, gmT2] = get_relaxation_times(3,'gm');
        [bloodT1, bloodT2] = get_relaxation_times(3,'blood');
        csfT1 = 3000;
        csfT2 = 300;
        T1(:,1) = wmT1;
        T1(:,2) = gmT1;
        T1(:,3) = bloodT1;
        T1(:,4) = csfT1;
        
        T2(:,1) = wmT2;
        T2(:,2) = gmT2;
        T2(:,3) = bloodT2;
        T2(:,4) = csfT2;
end
%
%
%% B1 distribution
B1 = fspecial('gaussian', dim, dim/6);
B1 = B1./(max(max(B1)));
%inver the distribution, to make it more like an expected B1 efficiency distribution
B1 = (-B1+1)/5;
B1 = (B1 + (1-max(B1(:))));
%% generate offset list
nPts = 48;
% add hard-coded sequence timings
minTR = 130; % ms, from sequence protocol card.
minTE = 32; % ms, from sequence protocol card.
%   offsets = generate_offset_list(lower,upper,nPts,seed,seq,saveFlag);
offsets = generate_offset_list([0 0 75 165],[50 50 105 195],nPts,[1 2 3 4],'SSEPI','noSave');

%% make dictionary
% Specify the list of dictionary parameters
T1s = 0:5:300;
T2s = 0:5:300;
%FAdevs = 0.75:0.01:1.25;
FAdevs = 1;
dictionaryParams(1,1:numel(T1s)) = T1s; % T1
dictionaryParams(2,1:numel(T2s)) = T2s ; % T2
dictionaryParams(3,1:numel(FAdevs)) = FAdevs ; % B1 fraction

nRuns = 1;
df = 0;
[dict, paramIndices] = compile_SE_dictionary_Bernstein(offsets, minTR, minTE, nRuns, df, dictionaryParams);

%% simulate acquired data
% [testFPtcMap, testT1map, testT2map] = generate_test_FP_data(dim1, dim2,nTimePts,offsetList,freqOffset,mapType,SNR)
SNR = 0;
simsig = generate_test_FP_data(T1,T2, minTR, minTE, nRuns,offsets,df,SNR,B1);


%% Compare dictionary with simulated signal
addpath '~/Documents/MATLAB/Ma2013'


SNRs = 10:10:100;
seqLs = 8:10:48;
for seqL = seqLs;
for SNR = SNRs;
    SNR
    for seed = 1:100;
        for n=1:size(simsig,2);
            % add random noise drawn from specific seed
            rng(seed);
            sig(:,n) = awgn(abs(simsig(:,n)),SNR);
        end
        %% Match Signal
        [matchout]=templatematch(dict(:,1:seqL),sig(1:seqL,:));
        matchedT1 = paramIndices(1,matchout);
        matchedT2 = paramIndices(2,matchout);
        matchedB1 = paramIndices(3,matchout);
        
        %% Calculate RMSE
        for nVoxel = 1:numel(T1)
            T1err(nVoxel,SNR,seqL,seed) = matchedT1(nVoxel) - T1(nVoxel);   % Errors
            T2err(nVoxel,SNR,seqL,seed) = matchedT2(nVoxel) - T2(nVoxel);   % Errors
        end
        
    end
    
end
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
cmin = 200;
cmax = 300;
%
subplot 341
imagesc(T1);
caxis([cmin cmax])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Synthetic T1'
%
subplot 345
imagesc(T2);
caxis([cmin cmax])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Synthetic T2'
%
subplot 342
imagesc(matchedT1);
caxis([cmin cmax])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Matched T1'
%
subplot 346
imagesc(matchedT2);
caxis([cmin cmax])
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
