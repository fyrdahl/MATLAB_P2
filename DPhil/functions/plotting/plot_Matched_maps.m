function plot_Matched_maps(matchedT1, matchedT2, matchedFAdev, M0fit_grad,varargin)


limits = varargin{1};
figure,
subplot 221
imagesc(matchedT1)
title 'Matched T1'
c = colorbar;
ylabel(c,'[ms]')
set(gca, 'FontSize',20)
axis square
colormap hot
caxis([limits(1,1), limits(1,2)])


subplot 222
imagesc(matchedT2)
title 'Matched T2'
c = colorbar;
ylabel(c,'[ms]')
set(gca, 'FontSize',20)
axis square
colormap hot
caxis([limits(2,1), limits(2,2)])

subplot 223
imagesc(matchedFAdev)
title 'Flip Angle Deviation'
c = colorbar;
ylabel(c,'Efficiency')
set(gca, 'FontSize',20)
axis square
colormap hot
caxis([limits(3,1), limits(3,2)])

subplot 224
imagesc(M0fit_grad)
title('Modulated Proton Density')
c = colorbar;
ylabel(c,['M_0','R'])
set(gca, 'FontSize',20)
axis square
colormap hot


end