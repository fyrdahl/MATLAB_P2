%% Script to Simulate IR-bSSFP data and assign properties
clear all
%
%% Generate offset list
%
disp 'Generating Sequence Parameter List...'
nPts = 48;

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

%% Make synthetic data
%
disp 'Making Synthetic Data...'

SNR = 75;
addpath('~/Documents/MATLAB/Ma2013')
dim = 64;
T1 = 282;
T2 = 214;
if dim > 1;
    T1 = phantom('Modified Shepp-Logan',dim)*T1;
    T2 = phantom('Modified Shepp-Logan',dim)*T2;
    T1 = T1(:);
    T2 = T2(:);
end
T1 = awgn(T1,10);
T2 = awgn(T2,10);
T1(T1<1) =0;
T2(T2<1) =0;
%
for n=1:numel(T1);
    simsig(:,n)=makeMRFdictionary(RFpulses ,offsets(:,1) ,T1(n), T2(n), df);
    signoise(:,n) = awgn(simsig(:,n),SNR);
end
sig = abs(signoise);
disp 'Making Synthetic Data... Done.'

%
%% Create Dictionary
disp 'Creating Signal Dictionary...'

dictT1 = [0:5:300];
dictT2 = [0:5:300];
index = 1;
for n = 1:length(dictT1)
    for m = 1:length(dictT2)
        tempdict=makeMRFdictionary(RFpulses ,offsets(:,1) ,dictT1(n), dictT2(m), df);
        dict(index,:) = abs(tempdict);
        lookUpTable(index,1)=dictT1(n);
        lookUpTable(index,2)=dictT2(m);
        index = index + 1;
    end
end

disp 'Creating Signal Dictionary... Done.'
%% Match Signal
disp 'Finding Best Matches...'
[matchout]=templatematch(dict,sig);
matchedT1 = lookUpTable(matchout,1);
matchedT2 = lookUpTable(matchout,2);

% RMSE 
for nVoxel = 1:numel(matchedT1)
T1err(nVoxel) = matchedT1(nVoxel) - T1(nVoxel);   % Errors
T2err(nVoxel) = matchedT2(nVoxel) - T2(nVoxel);   % Errors
end
T1RMSE = sqrt(mean((T1err).^2));  % Root Mean Squared Error
T2RMSE = sqrt(mean((T2err).^2));  % Root Mean Squared Error



disp 'Finding Best Matches... Done.'
%% Plot Results
disp 'Plotting Results...'
T1 = reshape(T1,dim,dim);
T2 = reshape(T2,dim,dim);
matchedT1 = reshape(matchedT1,dim,dim);
matchedT2 = reshape(matchedT2,dim,dim);
%
%
figure('name',['IRbSSFP - with Gaussian Noise (SNR = ',num2str(SNR),' nPts = ',num2str(nPts),')'])
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
c=colorbar;
colormap hot
ylabel(c,'Error [ms]')
title(['T1 Error (RMSE = ',num2str(T1RMSE),')'])
axis square

%
subplot 247
imagesc(abs(matchedT2-T2))
c=colorbar;
colormap hot
ylabel(c,'Error [ms]')
title(['T2 Error (RMSE = ',num2str(T2RMSE),')']) 
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


disp 'Plotting Results... Done.'