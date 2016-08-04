function plot_best_match_versus_data(bestMatch, data, data_y, data_x, M0fit_grad)
%% plot_best_match_versus_data(bestMatch, data, data_y, data_x, M0fit_grad)
% A function to plot the acquired signal against the simulated signal.
% This function can be used to check that the M0 measurements (obtained from the gradient of this plot)
% are correct.
%
% Receives 3D matrices.
%%
if data_x <= 0 || data_y <=0
    error('Coordinates must be greater than zero')
end
figure('name','Acquired Data Timecourse Versus Best Matched Entry')
plot(squeeze(bestMatch(data_x,data_y,:)),squeeze(data(data_x,data_y,:)),'x')
xlabel 'Best Match [a.u.]'
ylabel 'Acquired Data [a.u.]'
hold on
x = 0.01:0.001:max(squeeze(bestMatch(data_x,data_y,:)));
plot(x,x*squeeze(M0fit_grad(data_x,data_y)),'r')

%nPts = size(data,3);
%plot
end