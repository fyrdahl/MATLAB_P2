function plot_simulated_signal(M, Mxy, imageTimes)
% PLOT_SIMULATED_SIGNAL(M, MXY, IMAGETIMES)
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
%
% Function to plot magnetisation over time.
%
%   M is a 3-by-n dimensional matrix, containing the magnetisation
%   evolution in x, y and z.
%
%   Mxy is a m-element vector, where m is the number of images acquired.
%
%   IMAGE_TIMES is a m-element vector, where m is the number of images
%   acquired. It stores the times at which the images were acquired.

figure('name',['Simulated Magnetisation. ',datestr(now,30)])
hold on
plot(M(1,:))
plot(M(2,:))
plot(M(3,:))
plot(abs(complex(M(1,:),M(2,:))),'--k')
for i = 1:size(Mxy,1)
plot(imageTimes(i,:), Mxy(i,:),'+')
end
xlabel 'Time [ms]'
ylabel 'Magnetisation [a.u.]'
legend('Mx', 'My', 'Mz' ,'Mtransverse','Image Signal','location','best')
set(gca,'fontsize',20)

end