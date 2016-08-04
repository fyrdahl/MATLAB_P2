function plot_FP_parameters(fingerprintLists)
<<<<<<< HEAD
%% PLOT_FP_PARAMETERS(FINGERPRINTLISTS)
% Author: <jack.allen@jesus.ox.ac.uk>
% 
% FINGERPRINTLISTS is a n-by-m matrix, where n is the number of
% measurements and m is the number of parameters.
%
% Plots the values of each column.
=======
%% Plot_FP_parameters
% jack.allen@jesus.ox.ac.uk
% University of Oxford
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1

list = squeeze(fingerprintLists);

figure
subplot 221
plot(list(:,1),'-*')
<<<<<<< HEAD
set(gca,'FontSize',24)
xlabel ('Image Index')
ylabel ('Repetition Time (TR) [ms]')
subplot 222
plot(list(:,2),'-*')
set(gca,'FontSize',24)
xlabel ('Image Index')
ylabel ('Echo Time (TE) [ms]')
subplot 223
plot(list(:,3),'-*')
set(gca,'FontSize',24)
xlabel ('Image Index')
ylabel ('Flip Angle 1 [degrees]')
subplot 224
plot(list(:,4),'-*')
set(gca,'FontSize',24)
xlabel ('Image Index')
ylabel ('Flip Angle 2 [degrees]')

=======
xlabel 'Index'
ylabel 'Repetition Time (TR) [ms]'
subplot 222
plot(list(:,2),'-*')
xlabel 'Index'
ylabel 'Echo Time (TE) [ms]'
subplot 223
plot(list(:,3),'-*')
xlabel 'Index'
ylabel 'Flip Angle 1 [degrees]'
subplot 224
plot(list(:,4),'-*')
xlabel 'Index'
ylabel 'Flip Angle 2 [degrees]'
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end