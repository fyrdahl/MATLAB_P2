%% Script to Simulate IR-bSSFP data and assign properties
clear all
%

plotFlag = 'noPlot';
%% Generate offset list
%
disp 'Generating Sequence Parameter List...'
nPts = 1000;
%
offsets = generate_offset_list([10 10],[15 15],nPts,[1 2],'SSFP','noSave');
%
% set RF pulses
% rotation determined by abs()
% phase determined by angle()
for n = 1:nPts
    if mod(n,2)==1 %odd time point
        RFreal(n) = offsets(n,2);
        RFimag(n) = 0;
    end
    if mod(n,2)==0 %even
        RFreal(n) = 0; %degrees
        RFimag(n) = offsets(n,2);
    end
end
RFreal = 2*pi*RFreal/360;
RFimag = 2*pi*RFimag/360;
RFpulses = complex(RFreal,RFimag);
%
% Set Frequency offset
df = 0;
disp 'Generating Sequence Parameter List... Done.'

%% Make Synthetic Phantom
SNR = 75;
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
        B1 = fspecial('gaussian', dim, dim/12);
        B1 = B1./(max(max(B1)));
        %inver the distribution, to make it more like an expected B1 efficiency distribution
        B1 = (-B1+1)/4;
        B1 = (B1 + (1-max(B1(:))));
        B1 = B1(size(B1,1)/2,:);
        
end
%

%% Create Dictionary
disp 'Creating Signal Dictionary...'
%
switch phantom
    case 'simpleBrain'
        dictT1 = [700:5:900, 1300:5:1700, 2900:5:3100];
        dictT2 = [50:5:90, 130:5:350] ;
end
index = 1;
% make dictionary using combinations of T1 and T2
for n = 1:length(dictT1)
    for m = 1:length(dictT2)
        if dictT1(n)>dictT2(m)
            tempdict=makeMRFdictionary(RFpulses ,offsets(:,1) ,dictT1(n), dictT2(m), df);
            dict(index,:) = abs(tempdict);
            lookUpTable(index,1)=dictT1(n);
            lookUpTable(index,2)=dictT2(m);
        end
        index = index + 1;
    end
end
%
disp 'Creating Signal Dictionary... Done.'

%
%% Simulate Acquired Data
for n=1:numel(T1map);
    for nB1 = 1:numel(B1)
        RFpulses
        tmpRFpulses = RFpulses*B1(nB1)
        pause
        simsig(:,n,nB1)=makeMRFdictionary(tmpRFpulses ,offsets(:,1) ,T1map(n), T2map(n), df);
    end
end
disp 'Making Synthetic Data... Done.'


%% iterate through different SNR and noise seeds
SNRs = 2;
seqLs = [100:100:1000];
for B1 = 1%:numel(B1)
for seqL = seqLs;
    for SNR = SNRs;
        SNR
        for seed = 1:100;
            for n=1:size(simsig,2);
                % add random noise drawn from specific seed
                rng(seed);
                sig(:,n,nB1) = awgn(abs(squeeze(simsig(:,n,nB1))),SNR);
            end
            %% Match Signal
            [matchout]=templatematch(dict(:,1:seqL),sig(1:seqL,:));
            matchedT1 = lookUpTable(matchout,1);
            matchedT2 = lookUpTable(matchout,2);
            
            %% Calculate RMSE
            for nVoxel = 1:numel(T1map)
                T1err(nVoxel,SNR,seqL,nB1,seed) = matchedT1(nVoxel) - T1map(nVoxel);   % Errors
                T2err(nVoxel,SNR,seqL,nB1,seed) = matchedT2(nVoxel) - T2map(nVoxel);   % Errors
            end
            
        end
        
    end
end
end
%%
T1RMSE = sqrt(mean((T1err).^2));  % Root Mean Squared Error
T2RMSE = sqrt(mean((T2err).^2));  % Root Mean Squared Error
%%
nSNR = 10;
nSeqL = 100;

%%
for nVoxel = 1:numel(T1map)
    for nSNR = 10:10:100
        for nSeqL = 100:100:1000
            T1stds(nVoxel,nSNR,nSeqL) = std((squeeze(T1err(nVoxel,nSNR,nSeqL,:))));
            T2stds(nVoxel,nSNR,nSeqL) = std((squeeze(T2err(nVoxel,nSNR,nSeqL,:))));
            T1means(nVoxel,nSNR,nSeqL) = mean((squeeze(T1err(nVoxel,nSNR,nSeqL,:))));
            T2means(nVoxel,nSNR,nSeqL) = mean((squeeze(T2err(nVoxel,nSNR,nSeqL,:))));
        end
    end
end
%% Plot Results
switch plotFlag
    case 'plot'
        disp 'Plotting Results...'
        matchedT1 = reshape(matchedT1,dim,dim);
        matchedT2 = reshape(matchedT2,dim,dim);
        
        figure
        x = repmat(SNR,numel(SNR),1);
        y = repmat(seqL,numel(seqL),1);
        z = squeeze(T1means(1,10:10:100,100:100:1000));
        e = squeeze(T1stds(1,10:10:100,100:100:1000));
        plot3d_errorbars(x,y',z,e)
        %        figure('name',['IRbSSFP Monte Carlo Results: T1'])
        %         for nVoxel = 1:numel(T1)
        %         errorbar(SNRs,T1means(nVoxel,SNRs),T1stds(nVoxel,SNRs),'-o')
        %         hold on
        %         end
        %
        %
        figure('name',['IRbSSFP - Gaussian Noise and B1 (SNR = ',num2str(SNR),' nPts = ',num2str(nPts),')'])
        %
        cmin = 200;
        cmax = 300;
        %
        subplot 341
        imagesc(T1map);
        caxis([cmin cmax])
        c=colorbar;
        colormap hot
        ylabel(c,'T1 [ms]')
        title 'Synthetic T1'
        %
        subplot 345
        imagesc(T2map);
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
        imagesc(abs(matchedT1-T1map))
        %caxis([0 300])
        c=colorbar;
        colormap hot
        ylabel(c,'T1 [ms]')
        title(['Error (RMSE = ',num2str(T1RMSE),')'])
        %
        subplot 347
        imagesc(abs(matchedT2-T2map))
        %caxis([0 300])
        c=colorbar;
        colormap hot
        ylabel(c,'T2 [ms]')
        title(['Error (RMSE = ',num2str(T2RMSE),')'])
        %
        subplot 344
        plot(T1map(:),matchedT1(:),'o')
        hold on
        plot(0:300,0:300)
        xlim([0 300])
        ylim([0 300])
        ylabel(c,'T1 [ms]')
        ylabel 'Matched T1'
        xlabel 'Synthetic T1'
        %
        subplot 348
        plot(T2map(:),matchedT2(:),'o')
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
    case 'noPlot'
        
end