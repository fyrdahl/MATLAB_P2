clear all
%

plotFlag = 'noPlot';
%% Generate offset list
%
nPts = 1000;
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
for nVoxel = 901:nPts
    rng(k);
    FAs(nVoxel) = 0.5*(10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
%
% TR lengths
TR = perlin_noise(nPts, 10,15,1.2);
figure, plot(TR)
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
RFreal = 2*pi*RFreal/360;
RFimag = 2*pi*RFimag/360;
RFpulses = complex(RFreal,RFimag);
figure, plot(abs(RFpulses))
%
% Set Frequency offset
df = 0;
disp 'Generating Sequence Parameter List... Done.'

%% Make Synthetic Phantom
addpath('~/Documents/MATLAB/Ma2013')
dim = 32;
refT1 = 282;
refT2 = 214;

phantomName = 'simpleBrain';

switch phantomName
    
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
%

%% Create Dictionary
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
end
%calculate size of dictionary
numEntries = 1;
for nT1 = 1:length(dictT1)
    for nT2 = 1:length(dictT2)
        if dictT1(nT1)>dictT2(nT2)
            for ndf = 1:length(dictdf)
                numEntries = numEntries + 1;
            end
        end
    end
end
% make dictionary using combinations of T1 and T2
index = 1;
simSpace = [length(dictT1), length(dictT2)];
numSims = prod(simSpace);
for idx = 1:numSims
    [nT1, nT2, ndf] = ind2sub(simSpace,idx);
   
        if dictT1(nT1)>dictT2(nT2)

                tempdict=makeMRFdictionary(RFpulses ,TR ,dictT1(nT1), dictT2(nT2), dictdf(ndf));
                IRbSSFPdict(index,:) = tempdict;
                IRbSSFPlookUpTable(index,1)=dictT1(nT1);
                IRbSSFPlookUpTable(index,2)=dictT2(nT2);
                IRbSSFPlookUpTable(index,3)=dictdf(ndf);
                index = index + 1;
        end

end
%
disp 'Creating Signal Dictionary... Done.'
toc
save('~/Documents/DPhil/MAT-files/IRbSSFPlut','IRbSSFPlookUpTablelookUpTable')



%% iterate through different SNR and noise seeds
numSeeds = 100;
B1 = 1;
SNRs = 100;
seqLs = [50,100:100:1000];
simSpace = [numel(T1map),numel(SNRs),numel(seqLs),numel(B1),numSeeds];
numSims = prod(simSpace);
IRbSSFPT1err = zeros(1,numSims);
IRbSSFPT2err = zeros(1,numSims);
matchout = zeros(1,numSims);
IRbSSFPT1matched = zeros(1,numSims);
IRbSSFPT2matched = zeros(1,numSims);
disp 'Calculate Errors...'
tic


T = 0;
for idx=1:numSims
    idx

    % extract parameter indices from general index
    [nVoxel, nSNR, nSeqLs, nB1, nSeed] = ind2sub(simSpace,idx);
     
    % make sythetic data
     % simulate flip angle efficiency
     tmpRF = RFpulses*B1(nB1);
     %full length timecourse simulation
     tmpsimsig=makeMRFdictionary(tmpRF ,TR ,T1map(nVoxel), T2map(nVoxel), 0);    
     % add noise to simulated signal
     tmpsimsig = awgn(tmpsimsig',SNRs(nSNR),0,nSeed);
     tmpsimsig = tmpsimsig(1:seqLs(nSeqLs));
     tmpsimsig=tmpsimsig./norm(tmpsimsig); %JA: fixed bug in original code (was not normalised in correct direction)

 
     % shorten dictionary 
     tmpdict = IRbSSFPdict(:,1:seqLs(nSeqLs));
     tmpdict=tmpdict./repmat(sqrt(sum(tmpdict.^2,2)),1,seqLs(nSeqLs));

    % Calculate inner product
    innerproduct=tmpdict*tmpsimsig;
    % Take the maximum value and return the index
    [~,matchout]=max(abs(innerproduct));

     %calculate error
     IRbSSFPT1err(idx)= abs(IRbSSFPlookUpTable(matchout,1)-T1map(nVoxel));
     IRbSSFPT2err(idx)= abs(IRbSSFPlookUpTable(matchout,2)-T2map(nVoxel));
     IRbSSFPT1matched(idx) = IRbSSFPlookUpTable(matchout,1);
     IRbSSFPT2matched(idx) = IRbSSFPlookUpTable(matchout,2);
end