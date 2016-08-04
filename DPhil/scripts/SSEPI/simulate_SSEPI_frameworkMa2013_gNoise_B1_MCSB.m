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

phantomName = 'simpleBrain';

switch phantomName
    
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
        
        
        %% Generate Synthetic B1 Distribution
        B1 = fspecial('gaussian', dim, dim/12);
        B1 = B1./(max(max(B1)));
        %inver the distribution, to make it more like an expected B1 efficiency distribution
        B1 = (-B1+1)/5;
        B1 = (B1 + (1-max(B1(:))));
        B1 = B1(size(B1,1)/2,:);
end
%
%

%% generate offset list
nPts = 96;
% add hard-coded sequence timings
minTR = 130; % ms, from sequence protocol card.
minTE = 32; % ms, from sequence protocol card.
%

% First cycle
for nVoxel = 1:250
    rng(nVoxel);
    FAs(nVoxel) = 10 + abs(sin(nVoxel*(2*pi)/500)*(50+(10*rand)));
end
% gap
for nVoxel = 251:310
    FAs(nVoxel) = 0;
end
% second cycle
k = 1;
for nVoxel = 311:560
    rng(k);
    FAs(nVoxel) = 0.5*(10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
% gap
for nVoxel = 561:590
    FAs(nVoxel) = 0;
end
% third cycle
k = 1;
for nVoxel = 591:840
    rng(k);
    FAs(nVoxel) = (10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
% gap
for nVoxel = 841:900
    FAs(nVoxel) = 0;
end
% fourth cycle
k = 1;
for nVoxel = 901:1000
    rng(k);
    FAs(nVoxel) = 0.5*(10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
%
%
% set RF pulses
% rotation determined by abs()
% phase determined by angle()
for nVoxel = 1:nPts
    if mod(nVoxel,2)==1 %odd time point
        RFreal(nVoxel) = FAs(nVoxel);
        RFimag(nVoxel) = 0;
    end
    if mod(nVoxel,2)==0 %even
        RFreal(nVoxel) = 0; %degrees
        RFimag(nVoxel) = FAs(nVoxel);
    end
    
end
pn = perlin_noise(1000, 0,5,1.2);
offsets(:,1) = pn(1:nPts);
offsets(:,2) = pn(1:nPts);
offsets(:,3) = abs(complex(RFreal,RFimag)); %just the rotation angle (degrees), no phase
offsets(:,4) = 2*abs(complex(RFreal,RFimag)); %just the rotation angle (degrees), no phase
offsets = generate_offset_list([0 0 75 165],[50 50 105 195],nPts,[1 2 3 4],'SSEPI','noSave')
figure, plot(offsets(:,1))
figure, plot(offsets(:,2))
figure, plot(offsets(:,3))
figure, plot(offsets(:,4))
%
% Set Frequency offset
df = 0;
disp 'Generating Sequence Parameter List... Done.'
%%  Make Synthetic Phantom
addpath('~/Documents/MATLAB/Ma2013')
dim = 32;
refT1 = 282;
refT2 = 214;

phantom = 'simpleBrain';

switch phantom
    
    case 'SL'
        if dim > 1;
            T1map = round(phantom('Modified Shepp-Logan',dim)*10);
            T2map = round(phantom('Modified Shepp-Logan',dim)*10);
            T1map(T1map<2) = 0;
            T2map(T2map<2) = 0;
            
            T1map(T1map==2) = 1*refT1;
            T1map(T1map==3) = 0.9*refT1;
            T1map(T1map==4) = 0.95*refT1;
            T1map(T1map==10) = 0.9*refT1;
            
            T2map(T2map==2) = 1*refT2;
            T2map(T2map==3) = 1.05*refT2;
            T2map(T2map==4) = 1*refT2;
            T2map(T2map==10) = 1.1*refT2;
            
            T1map = awgn(T1map,5);
            T2map = awgn(T2map,5);
            
            T1map(T1map<10) = 0;
            T2map(T2map<10) = 0;
            
            % Generate Synthetic B1 Distribution
            B1 = fspecial('gaussian', dim, dim/6);
            B1 = B1./(max(max(B1)));
            %inver the distribution, to make it more like an expected B1 efficiency distribution
            B1 = (-B1+1)/5;
            B1 = (B1 + (1-max(B1(:))));
            
        end
        
    case 'Brain' %Brain Values
        T1map = ones(dim);
        T2map = ones(dim);
        
        [wmT1, wmT2] = get_relaxation_times(3,'wm');
        [gmT1, gmT2] = get_relaxation_times(3,'gm');
        [bloodT1, bloodT2] = get_relaxation_times(3,'blood');
        csfT1 = 3000;
        csfT2 = 300;
        T1map(:,1:dim/4) = wmT1;
        T1map(:,dim/4:dim/2) = gmT1;
        T1map(:,dim/2:(3*dim/4)) = bloodT1;
        T1map(:,(3*dim/4):end) = csfT1;
        
        T2map(:,1:dim/4) = wmT2;
        T2map(:,dim/4:dim/2) = gmT2;
        T2map(:,dim/2:(3*dim/4)) = bloodT2;
        T2map(:,(3*dim/4):end) = csfT2;
        
        % Generate Synthetic B1 Distribution
        B1 = fspecial('gaussian', dim, dim/6);
        B1 = B1./(max(max(B1)));
        %inver the distribution, to make it more like an expected B1 efficiency distribution
        B1 = (-B1+1)/5;
        B1 = (B1 + (1-max(B1(:))));
        
    case 'simpleBrain' %Brain Values
        nTissue = 4;
        T1map = ones(1,nTissue);
        T2map = ones(1,nTissue);
        
        [wmT1, wmT2] = get_relaxation_times(3,'wm');
        [gmT1, gmT2] = get_relaxation_times(3,'gm');
        [bloodT1, bloodT2] = get_relaxation_times(3,'blood');
        csfT1 = 3000;
        csfT2 = 300;
        T1map(:,1) = wmT1;
        T1map(:,2) = gmT1;
        T1map(:,3) = bloodT1;
        T1map(:,4) = csfT1;
        
        T2map(:,1) = wmT2;
        T2map(:,2) = gmT2;
        T2map(:,3) = bloodT2;
        T2map(:,4) = csfT2;
        
        % Generate Synthetic B1 Distribution
%         B1 = fspecial('gaussian', dim/4, dim/16);
%         B1 = B1./(max(max(B1)));
%         %inver the distribution, to make it more like an expected B1 efficiency distribution
%         B1 = (-B1+1)/4;
%         B1 = (B1 + (1-max(B1(:))));
%         B1 = B1(size(B1,1)/2,:);
%         B1 = B1/max(B1)
%         
         B1 = 0.8:0.1:1.2;
end





%% make dictionary
% Specify the list of dictionary parameters

disp 'Creating Signal Dictionary...'
%
tic
switch phantomName
    case 'simpleBraindf'
        dictT1 = [100:20:2000, 2300:300:5000];% ms
        dictT2 = [20:5:100, 110:10:200, 400:200:3000] ;% ms
        dictdf = [-390:20:-270,...
            -250:10:-90,...
            -80:2:-42,...
            -40:40,...
            42:2:80,...
            90:10:250,...
            270:20:390]; % Hz
    case 'simpleBrain'
       % dictT1 = [500:5:2200, 2200:5:4000];
       % dictT2 = [0:5:500] ;
       % dictdf = 0; %Hz
        dictT1 = [700:5:900, 1100:5:1700, 2700:5:3300];
        dictT2 = [20:5:350] ;
        dictdf = 0;
        dictFAdevs = 1;
end


dictionaryParams(1,1:numel(dictT1)) = dictT1; % T1
dictionaryParams(2,1:numel(dictT2)) = dictT2 ; % T2
dictionaryParams(3,1:numel(dictFAdevs)) = dictFAdevs ; % B1 fraction
%
%calculate size of dictionary
numEntries = 1;
for nT1 = 1:numel(dictT1)
    for nT2 = 1:numel(dictT2)
        if dictT1(nT1)>=dictT2(nT2)
            for ndf = 1:length(dictdf)
                for nFAdev = 1:numel(dictFAdevs)
                numEntries = numEntries + 1;
                end
            end
        end
    end
end
nRuns = 1;
df = 0;
%[dict,pI1,pI2,pI3, lookUpTable] = compile_SE_dictionary_Bernstein(offsets,nPts, minTR, minTE, nRuns, df, dictionaryParams);
[dictionary, paramIndices] = compile_SE_dictionary_Bernstein(offsets, minTR, minTE, nRuns, df, dictionaryParams);
lookUpTable = paramIndices';
%% simulate acquired data
B1 = 1;
numSeeds = 100;
SNRs = 50;
seqLs = [10,30,70,96] %:10:100;
simSpace = [numel(T1map),numel(SNRs),numel(seqLs),numel(B1),numSeeds];
numSims = prod(simSpace);
% pre-allocate results matrices
SSEPIT1err = zeros(1,numSims);
SSEPIT2err = zeros(1,numSims);
matchout = zeros(1,numSims);
T1matched = zeros(1,numSims);
T2matched = zeros(1,numSims);
tmpoffsets = offsets;
disp 'Calculate Errors...'
tic

T = 0;

for idx=1:numSims
    idx

    % extract parameter indices from general index
    [nVoxel, nSNR, nSeqLs, nB1, nSeed] = ind2sub(simSpace,idx);
    
     % make sythetic data
     % simulate flip angle efficiency
     tmpoffsets(:,3) = offsets(:,3)*B1(nB1);
     tmpoffsets(:,4) = offsets(:,4)*B1(nB1);
     tmpsimsig = sim_SE_bernstein(T1map(nVoxel), T2map(nVoxel), minTR, minTE, tmpoffsets,nRuns,df);
     tmpsimsig = awgn(tmpsimsig,SNRs(nSNR),0,nSeed); 
     tmpsimsig = tmpsimsig(1:seqLs(nSeqLs));
     tmpsimsig = tmpsimsig./norm(tmpsimsig);

     %normalise dictionary
     tmpdictionary = dictionary(:,1:seqLs(nSeqLs));
     tmpdictionary=tmpdictionary./repmat(sqrt(sum(tmpdictionary.^2,2)),1,seqLs(nSeqLs));
     
    % Calculate inner product
    innerproduct=tmpdictionary*tmpsimsig;
    % Take the maximum value and return the index
    [~,matchout]=max(abs(innerproduct));

     %calculate error
     SSEPIT1err(idx)= abs(lookUpTable(matchout,1)-T1map(nVoxel));
     SSEPIT2err(idx)= abs(lookUpTable(matchout,2)-T2map(nVoxel));
     % store matched T1 and T2
     T1matched(idx) = lookUpTable(matchout,1);
     T2matched(idx) = lookUpTable(matchout,2);
end
T1matched = reshape(T1matched,simSpace);
T2matched = reshape(T2matched,simSpace);
SSEPIT1err = reshape(SSEPIT1err,simSpace);
SSEPIT2err = reshape(SSEPIT2err,simSpace);

disp 'Calculate Errors... Done.'
toc
%
%%
for nB1 = 1:numel(B1)
for nVoxel = 1:numel(T1)
    for nSNR = 1:numel(SNRs)
        for nSeqL = 1:numel(seqLs)
            T1stds(nVoxel,nSNR,nSeqL,nB1) = std((squeeze(SSEPIT1err(nVoxel,nSNR,nSeqL,nB1,:))));
            T2stds(nVoxel,nSNR,nSeqL,nB1) = std((squeeze(SSEPIT2err(nVoxel,nSNR,nSeqL,nB1,:))));
            T1means(nVoxel,nSNR,nSeqL,nB1) = mean((squeeze(SSEPIT1err(nVoxel,nSNR,nSeqL,nB1,:))));
            T2means(nVoxel,nSNR,nSeqL,nB1) = mean((squeeze(SSEPIT2err(nVoxel,nSNR,nSeqL,nB1,:))));
        end
    end
end
end

%% Plot Results
figure
x = repmat(SNRs,numel(seqLs),1);
y = repmat(seqLs,numel(SNRs),1);
z = squeeze(T1means(1,SNRs,seqLs));
e = squeeze(T1stds(1,SNRs,seqLs));
plot3d_errorbars(x,y',z,e,'surf')






%%
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
