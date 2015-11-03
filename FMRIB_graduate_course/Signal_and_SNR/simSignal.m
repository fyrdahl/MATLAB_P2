% written by Daniel Gallichan 05/2007
% modified by Karla Miller 05/2008

% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

if ~exist('T1','var') | ~exist('TR','var') | ~exist('flipAngle','var') | ~exist('tMax','var')
    disp('ERROR: T1, TR, flipAngle and tMax must all be defined')
    return
end
    
flip = flipAngle;

t = linspace(0,tMax,5000);
dt = t(2)-t(1);

if TR > tMax
    disp('ERROR: End Time must be longer than TR')
    return
end

M = zeros(size(t));
M(1) = 1; %cos(flip*pi/180);
tFlips = TR:TR:tMax;
iFlips = round(interp1(t,1:length(t),tFlips));
for iT = 2:length(t)
    M(iT) = M(iT-1) + ( (1-M(iT-1))/T1 )*dt;
    if any(iT == iFlips)
        M(iT) = M(iT)*cos(flip*pi/180);
    end
end

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])

signal = zeros(size(iFlips));
subplot(212);
for iF = 1:length(iFlips)
    signal(iF) = M(iFlips(iF)-1)*sin(flip*pi/180);
    hleg(2) = line([t(iFlips(iF)) t(iFlips(iF))],[0 signal(iF)],'color','r','linewidth',8);
end
axis([0 1.05*tMax 0 1.05]); grid on
xlabel('Time (ms)'); ylabel('M_{xy} (signal)'); 


subplot(211);
hleg(1) = plot(t,M,'linewidth',2);
axis([0 1.05*tMax 0 1.05]); grid on
ylabel('M_z'); 
title(['{\bf Tissue:} T1 = ' num2str(T1) 'ms   {\bf Scan:} TR = ' num2str(TR) 'ms, Flip Angle = ' num2str(flip) ' degrees        {\bf Resulting signal:} ' num2str(signal(end),3)])
fontScale(1.3)

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')

