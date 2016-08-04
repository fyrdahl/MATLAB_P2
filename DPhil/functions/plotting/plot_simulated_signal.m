function plot_simulated_signal(M, Mxy, imageTimes)
<<<<<<< HEAD
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
=======
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1

figure('name',['Simulated Magnetisation. ',datestr(now,30)])
hold on
plot(M(1,:))
plot(M(2,:))
plot(M(3,:))
plot(abs(complex(M(1,:),M(2,:))),'--k')
for i = 1:size(Mxy,1)
<<<<<<< HEAD
plot(imageTimes(i,:), Mxy(i,:),'+')
=======
plot(imageTimes, Mxy(i,:),'+')
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end
xlabel 'Time [ms]'
ylabel 'Magnetisation [a.u.]'
legend('Mx', 'My', 'Mz' ,'Mtransverse','Image Signal','location','best')
set(gca,'fontsize',20)

end