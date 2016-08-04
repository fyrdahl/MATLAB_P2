%% Script to Simulate IR-bSSFP data and assign properties
clear all
%
%% Generate offset list
%
nTimePts = 1000;

offsets = generate_offset_list([10 10],[15 15],nTimePts,[1 2],'SSFP','noSave');

%
% set RF pulses
% rotation determined by abs()
% phase determined by angle()
for n = 1:nTimePts
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

%% Make synthetic data
%
addpath('~/Documents/MATLAB/Ma2013')
dim = 10;
T1 = 282;
T2 = 214;
T1 = phantom('Modified Shepp-Logan',dim)*T1;
T2 = phantom('Modified Shepp-Logan',dim)*T2;
T1 = T1(:);
T2 = T2(:);
%
for n=1:numel(T1);
sig(:,n)=makeMRFdictionary(RFpulses ,offsets(:,1) ,T1(n), T2(n), df);
end
sig = abs(sig);
%
%% Create Dictionary
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
%% Match Signal
[matchout]=templatematch(dict,sig);
matchedT1 = lookUpTable(matchout,1);
matchedT2 = lookUpTable(matchout,2);
%% Plot Results
T1 = reshape(T1,dim,dim);
T2 = reshape(T2,dim,dim);
matchedT1 = reshape(matchedT1,dim,dim);
matchedT2 = reshape(matchedT2,dim,dim);
%
%
figure('name','IRbSSFP MRF Framework Simulation')
subplot 221
imagesc(matchedT1);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Matched T1'
%
subplot 222
imagesc(matchedT2);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Matched T2'
%
subplot 223
imagesc(T1);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Synthetic T1'
%
subplot 224
imagesc(T2);
caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Synthetic T2'
%
figure
subplot 221
imagesc(matchedT1-T1)
%caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T1 [ms]')
title 'Matched T1 - Synthetic T1'
%
subplot 222
imagesc(matchedT2-T2)
%caxis([0 300])
c=colorbar;
colormap hot
ylabel(c,'T2 [ms]')
title 'Matched T2 - Synthetic T2'
%
subplot 223
plot(T1(:),matchedT1(:),'o')
hold on
plot(0:300,0:300)
xlim([0 300])
ylim([0 300])
ylabel(c,'T1 [ms]')
ylabel 'Matched T1'
xlabel 'Synthetic T1'
%
subplot 224
plot(T2(:),matchedT2(:),'o')
hold on
plot(0:300,0:300)
xlim([0 300])
ylim([0 300])
ylabel(c,'T2 [ms]')
ylabel 'Matched T2'
xlabel 'Synthetic T2'








