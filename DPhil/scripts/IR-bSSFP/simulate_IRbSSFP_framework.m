%% Script to Simulate IR-bSSFP data and assign properties
clear all
%
%% Generate offset list
%
nTimePts = 1000;
% First cycle
for n = 1:250
    rng(n);
    FAs(n) = 10 + abs(sin(n*(2*pi)/500)*(50+(10*rand)));
end
% gap
for n = 251:310
    FAs(n) = 0;
end
% second cycle
k = 1;
for n = 311:560
    rng(k);
    FAs(n) = 0.5*(10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
% gap
for n = 561:590
    FAs(n) = 0;
end
% third cycle
k = 1;
for n = 591:840
    rng(k);
    FAs(n) = (10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
% gap
for n = 841:900
    FAs(n) = 0;
end
% fourth cycle
k = 1;
for n = 901:nTimePts
    rng(k);
    FAs(n) = 0.5*(10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
%
% TR lengths
TR = perlin_noise(nTimePts, 10,15,3);
%
% set RF pulses
% rotation determined by abs()
% phase determined by angle()
for n = 1:nTimePts
    if mod(n,2)==1 %odd time point
        RFreal(n) = FAs(n);
        RFimag(n) = 0;
    end
    if mod(n,2)==0 %even
        RFreal(n) = 0; %degrees
        RFimag(n) = FAs(n);
    end
    
end
RFreal = 2*pi*RFreal/360;
RFimag = 2*pi*RFimag/360;
RFpulses = complex(RFreal,RFimag);
%RFpulses(nTimePts/2-100:nTimePts/2) = 0;
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
sig(:,n)=makeMRFdictionary(RFpulses ,TR ,T1(n), T2(n), df);
end
sig = abs(sig);
%
%% Create Dictionary
dictT1 = [0:5:300];
dictT2 = [0:5:300];
index = 1;
for n = 1:length(dictT1)
    for m = 1:length(dictT2)
        tempdict=makeMRFdictionary(RFpulses ,TR ,dictT1(n), dictT2(m), df);
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

figure
subplot 221
imagesc(matchedT1);
subplot 222
imagesc(matchedT2);
subplot 223
imagesc(T1);
subplot 224
imagesc(T2);





