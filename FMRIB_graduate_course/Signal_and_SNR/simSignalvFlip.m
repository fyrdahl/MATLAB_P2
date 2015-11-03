if ~exist('T1','var') || ~exist('TR','var')
    disp('ERROR: T1 and  TR must be specified')
    return
end

% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

figure
theta = linspace(0,90,100);
sig = (1-exp(-TR/T1)).*sin(theta*pi/180)./(1-cos(theta*pi/180)*exp(-TR/T1));
plot(theta,sig,'linewidth',2); 
xlabel('Flip Angle (\circ)')
ylabel('Relative signal')
grid on
title(['{\bf Tissue:} T1 = ' num2str(T1) 'ms   {\bf Scan:} TR = ' num2str(TR) 'ms'])
fontScale(1.2)

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')
