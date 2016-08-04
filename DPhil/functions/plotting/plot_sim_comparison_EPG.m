function plot_sim_comparison_EPG(data, T1, T2, offsets)
% function to compare simulations with acquired data

% simulate the signal using EPGs
simsig = sim_SE_EPG(complex(deg2rad(offsets(:,3)),0),complex(deg2rad(offsets(:,4)),0),offsets(:,1)/1000,offsets(:,2)/1000,T1/1000,T2/1000)
% normalise 
simsig = abs(simsig);
simsig = simsig(:)./norm(simsig);

% get ROI
fig = figure;
imagesc(squeeze(data(:,:,1)));
rect = getrect(fig);


% Plot acquired signal
figure
for n = round(rect(1)):(round(rect(1))+round(rect(3)))
    for m = round(rect(2)):(round(rect(2))+round(rect(4)))
    y = squeeze(data(n,m,1:numel(simsig)));
    y = y(:)./norm(y);
    plot(y,'-o')
    hold on
   end
end

% plot simulated signal
plot(simsig,'b-*','LineWidth',3,'MarkerSize',15)

% figure settings
xlabel 'TE Index'
ylabel 'Normalised Signal'
set(gca, 'FontSize',18)

end