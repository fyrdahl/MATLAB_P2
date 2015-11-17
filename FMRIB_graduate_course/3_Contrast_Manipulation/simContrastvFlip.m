% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

if ~exist('T1pair','var') || ~exist('TR','var')
    disp('ERROR: T1pair and TR must be specified')
    return
end

if length(T1pair)~=2
    disp('ERROR: T1pair must contain 2 values')
    return
end

figure
set(gcf,'Position',[    50   679   740   409])
theta = linspace(0,90,1000);
sig = zeros(length(theta),2);
for iT1 = 1:2
    sig(:,iT1) = (1-exp(-TR/T1pair(iT1))).*sin(theta*pi/180)./(1-cos(theta*pi/180)*exp(-TR/T1pair(iT1)));
end

plot(theta,sig(:,1),'linewidth',2); hold all
plot(theta,sig(:,2),'linewidth',2); 
h1 = plot(theta,abs(sig(:,1)-sig(:,2)));
set(h1,'linewidth',2)
xlabel('Flip Angle (degrees)')
ylabel('Relative signal')
grid on
title(['TR = ' num2str(TR) 'ms'])
legend(['T1 = ' num2str(T1pair(1)) 'ms'], ['T1 = ' num2str(T1pair(2)) 'ms'], 'Contrast');
fontScale(1.4)

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')
