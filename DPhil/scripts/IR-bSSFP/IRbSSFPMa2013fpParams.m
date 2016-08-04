%% Script to Simulate IR-bSSFP data and assign properties
clear all
%
%% Generate offset list
%
nPts = 1000;
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
for n = 901:nPts
    rng(k);
    FAs(n) = 0.5*(10 + abs(sin(k*(2*pi)/500)*(50+(10*rand))));
    k = k+1;
end
%
% TR lengths
TR = perlin_noise(nPts, 10,15,3);
%
% set RF pulses
% rotation determined by abs()
% phase determined by angle()
for n = 1:nPts
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

offsetList(:,1) = TR;
offsetList(:,2) = rad2deg(abs(RFpulses));

fIDssfp = fopen(['~/Documents/MATLAB/DPhil/offsetLists/NEW_SSFP_Fingerprint_List.txt'],'a');
fprintf(fIDssfp,['LIST_OFFSET \r\n']);
fprintf(fIDssfp,[num2str(nPts),'\r\n']);
for i = 1:size(offsetList,1)
    fprintf(fIDssfp,[num2str(offsetList(i,1)),' ']);
    fprintf(fIDssfp,num2str(offsetList(i,2)));
    fprintf(fIDssfp,'\r\n');
end
