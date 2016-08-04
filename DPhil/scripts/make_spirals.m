%% Make Spirals

addpath(genpath('~/Documents/MATLAB/vdspiral'))

%% 1. Initialise

%	Examples using vds.m and vdsmex.m
%
%	This is a 16-interleave spiral that achieves
%	1mm resolution, and density decreases linearly
%	from supporting 24cm FOV at |k|=0 to 12cm FOV
%	at |k|=maximum.
%
%smax = 15000;	 % 150 T/m/s
smax = 18000;	 % 180 T/m/s

%
%((24/1000)*10000)/100
gmax = 2.4;	 % G/cm (system gmax = 24mT/m)

%T = .000004;	 % Seconds
T = .00001;	 % 10 microseconds
%
N = 48;		 % Interleaves
nInterleaves = N;

Fcoeff = [24 -12]; 	% FOV decreases linearly from 24 to 12cm.

res = 1;
rmax = 5/res;		% cm^(-1), corresponds to 1mm resolution.

%% 2. Generate initial spiral trajectory

disp('Calculating Gradient...');
[k,g,s,time,r,theta,FOV] = vds(smax,gmax,T,N,Fcoeff,rmax);

g = [ real(g(:)), imag(g(:))];

nSamplePts = size(g,1); %number of ADC columns needed

%% 3. Append gradients to return to (0,0) in k-space
count = size(g,1);
deltatime = time(count) -  time(count-1);
%
start = count;
while g(count,1) ~= 0
    g(count+1,1) = g(count,1) - smax*T;
    g(count+1,2) = g(count,2);
    if g(count,1) < (smax*T)
        g(count+1,1) =0;
    end
    time(count+1) = time(count) + T;
    
    count = count + 1;
    
    
end

while g(count,2) ~= 0
    g(count+1,2) = g(count,2) - smax*T;
    if g(count,2) < (smax*T)
        g(count+1,2) =0;
    end
    time(count+1) = time(count) + T;
    count = count + 1;
end

g(end+1,:) = 0;
time(end+1) = time(end) + T;
count = size(g,1);
g = complex(g(:,1),g(:,2));
gamma = 4258;
k = cumsum(g)*gamma*T;

%
nTpts = 100;
deltat = nTpts*T;
deltagx = (real(k(end))/gamma)/(deltat);
deltagy = (imag(k(end))/gamma)/(deltat);
n = 1;
residualgx(n) = 0;
residualgy(n) = 0;
%

for n = 1:(nTpts/2)
    residualgx(n+1) = residualgx(n) - deltagx/(deltat/(4*T));
    residualgy(n+1) = residualgy(n) + deltagy/(deltat/(4*T));
    
    n = n + 1;
end
for n= ((nTpts/2)+1):(nTpts)
    residualgx(n+1) = residualgx(n) + deltagx/(deltat/(4*T));
    residualgy(n+1) = residualgy(n) - deltagy/(deltat/(4*T));
    
    n = n + 1;
end

residualg = complex(residualgx,residualgy);
totalg = [g',residualg];
totalk = cumsum(totalg)*gamma*T;

r = abs(totalg);
theta = angle(totalg);

%% 4. Make additional interleaves
%
complexgradients = zeros(N,numel(r));
angleIncrement = (2*pi)/nInterleaves; % this is 7.5 degrees for N = 48
for spiralInd = 1:nInterleaves
    theta = theta+angleIncrement;
    complexgradients(spiralInd,:) = r.*exp(1i*(theta));
end
complexgradients(:,end+1) = 0; %gradients must end with a zero
%

%%  Save Gradient Waveforms as text file
%
spiralgradients = [real(complexgradients);imag(complexgradients)];
% normalise gradients for use in pulse sequence
spiralgradients = spiralgradients/max(max(spiralgradients));
%
% all pairs of gradient waveforms
fID = fopen('~/Documents/MATLAB/DPhil/spirals.txt','w'); %write a text file in the same folder as this script
for nRow = 1:(size(spiralgradients,1)-1)
    for nCol= 1:(size(spiralgradients,2)-1)
        fprintf(fID,'%.20f ',spiralgradients(nRow,nCol));
    end
    fprintf(fID,'%.20f',spiralgradients(nRow,(size(spiralgradients,2)))); % last entry in the spiral
    
    fprintf(fID,'\r\n'); %use /r/n to avoid windows formatting problems
end
%   last spiral
for nCol= 1:(size(spiralgradients,2)-1)
    fprintf(fID,'%.20f ',spiralgradients(size(spiralgradients,1),nCol));
end
fprintf(fID,'%.20f',spiralgradients(size(spiralgradients,1),(size(spiralgradients,2)))); % last entry in the spiral
%
%
% Save one pair of gradient waveforms to a text file as an example
fID2 = fopen('~/Documents/MATLAB/DPhil/singleSpiral.txt','w');
for nCol= 1:size(spiralgradients,2)
    fprintf(fID2,'%.20f ',spiralgradients(1,nCol));
end


%% calculate kspace trajectories from gradients
for spiralInd = 1:size(complexgradients,1)
    spiralkspace(spiralInd,:) = cumsum(complexgradients(spiralInd,:))*gamma*T;
end
%

%% Plot particular k-space trajectory
%
figure,
plot(spiralkspace(1,:))
%
figure,
plot(real(spiralkspace(1,:)))
hold on
plot(imag(spiralkspace(1,:)))
%
figure
for n = 1:48
    plot(spiralkspace(n,:),'.')
    hold on
end
%
%% Plot all k-space trajectories
%{
v = VideoWriter('spirals.avi');
if nargin ==2
    v.FrameRate = 2;
else
    v.FrameRate = 2;
end
open(v)


figure
subplot 121
plot(real(squeeze(complexspirals(1,:))), imag(squeeze(complexspirals(1,:))))
title 'Spiral 1'
set(gca,'FontSize',18)
axis square
xlim([-5 5])
ylim([-5 5])

subplot 122
title 'Overlay of All Spirals'
for spiralInd = 1:N
    hold on
    plot(real(squeeze(complexspirals(spiralInd,:))), imag(squeeze(complexspirals(spiralInd,:))),'-');
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'FontSize',18)
    axis square
    frame = getframe;
    writeVideo(v,frame);
end

close(v);

%}
%% Write to text file




