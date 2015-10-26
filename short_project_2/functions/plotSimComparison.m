function plotSimComparison(Mxy, data,compartmentCenters,compartmentCentersList, nTimeCoursePts,sliceNumber,phantomName,workingdir)

simCom = figure; hold on
% title (['fingerprint offset list ', num2str(offsetListNum)])
plot(Mxy(:)./Mxy(1),'x','MarkerSize',20)
y = zeros(3,nTimeCoursePts);

%% signal from each compartment
%title (['Phantom: ',phantomName,', Offset list: ',num2str(offsetListNum),', compartment center coords list: ',num2str(compartmentCentersList)]);
plotNumCompartments = 6;
for n = 1:plotNumCompartments
    %normStdBG = (std(background(:)))/y(n,1);
    data(n,:) = data(n,:)/data(n,1 );
    % residuals = y(n,:) - ySim;
    plot(data(n,:),'.') 
    % errorbar(y(n,:),repmat(normStdBG,1,24),'.' );  
end
ylim([0 1.1])
legend ({'Simulated Signal', 'Sample Pixel 1','Sample Pixel 2','Sample Pixel 3','Sample Pixel 4','Sample Pixel 5','Sample Pixel 6'},'Position',[0.6, 0.7, 0.1,0.1])
xlabel 'TE Index'
ylabel 'Normalised Signal'   
%% figure; plot(residuals,'+')

end