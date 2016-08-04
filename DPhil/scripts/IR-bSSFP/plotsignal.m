nPts = 256;
%
addpath(genpath('~/Documents/MATLAB/Ma2013'))
for k = 100
offsets = generate_offset_list([10 10],[k 10],nPts,[1 2],'SSFP','save',1);
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


nTissue = 4;
T1s = ones(1,nTissue);
T2s = ones(1,nTissue);

[wmT1, wmT2] = get_relaxation_times(3,'wm');
[gmT1, gmT2] = get_relaxation_times(3,'gm');
[bloodT1, bloodT2] = get_relaxation_times(3,'blood');
csfT1 = 3000;
csfT2 = 300;
T1s(:,1) = wmT1;
T1s(:,2) = gmT1;
T1s(:,3) = bloodT1;
T1s(:,4) = csfT1;

T2s(:,1) = wmT2;
T2s(:,2) = gmT2;
T2s(:,3) = bloodT2;
T2s(:,4) = csfT2;

figure
for n = 1:nTissue
    T1 = T1s(n);
    T2 = T2s(n);
    tempdict=makeMRFdictionary(RFpulses ,offsets(:,1) ,T1, T2, df);
    plot(abs(tempdict))
    hold on
end
 line([48 48],[0 1])
  line([2*48 2*48],[0 1])
    line([3*48 3*48],[0 1])
        line([4*48 4*48],[0 1])
            line([256 256],[0 1])
            
end