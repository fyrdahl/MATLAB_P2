<<<<<<< HEAD
function plot_sim_comparison(data, T1, T2, offsets,nRuns,df)
% plot_sim_comparison(Mxy,data,nPts,phantomName,savingdir,offsetListNum, sliceNumber)
%


%% simulate the signal using Bernstein mathematical description
%function [Msig, M] = sim_SE_bernstein(T1, T2, minTR, minTE, offsets,nRuns,df)
Msig = sim_SE_EPG(T1, T2, TRmin, TEmin, offsets,nRuns, df);
% normalise 
simSig = Msig/norm(Msig);


fig = figure;
imagesc(squeeze(data(:,:,1)));
rect = getrect(fig);

%% Plot simulated signal
figure

%% Plot acquired signal
nTCs = size(data,1);
nPts = size(data,3);
y = zeros(nTCs,nPts);
    
for n = round(rect(1)):(round(rect(1))+round(rect(3)))
    for m = round(rect(2)):(round(rect(2))+round(rect(4)))
    y = squeeze(data(n,m,1:numel(simSig)))
    y = y(:)/y(1);
    plot(y,'-o')
    hold on
   end
end

    % errorbar(y(n,:),repmat(normStdBG,1,24),'.' );
%end
plot(simSig,'b-*','LineWidth',3,'MarkerSize',15)
hold on

ylim([0 max(max(simSig))])
%legend ({'Simulated Signal', 'Sample Pixel 1','Sample Pixel 2','Sample Pixel 3','Sample Pixel 4','Sample Pixel 5','Sample Pixel 6'},'location','best')
xlabel 'TE Index'
ylabel 'Normalised Signal'

set(gca, 'FontSize',18)

%%
%savefig([savingdir,'/figures/compareSimwithData_Phantom_',phantomName,'__Offset_list_',num2str(offsetListNum),'_compartmentcentercoordslist:',num2str(compartmentCentersList),'.fig'])
%matlab2tikz('figurehandle',simCom,'filename',[savingdir,'/DTC_report/',phantomName,'simCom',num2str(offsetListNum),'slice',num2str(sliceNumber)],'height', '\figureheight', 'width', '\figurewidth')
% figure; plot(residuals,'+')
=======
function plotSimComparison(Mxy,data,nPts,phantomName,savingdir,compartmentCentersList,offsetListNum, sliceNumber)

simCom = figure; hold on
% title (['fingerprint offset list ', num2str(offsetListNum)])
plot(Mxy(:)./Mxy(1),'x','MarkerSize',20)
nTCs = size(data,1);
y = zeros(nTCs,nPts);

%% signal from each compartment
%title (['Phantom: ',phantomName,', Offset list: ',num2str(offsetListNum),', compartment center coords list: ',num2str(compartmentCentersList)]);

for n = 1:nTCs  
        y(n,:) = data(n,:);
    %normStdBG = (std(background(:)))/y(n,1);
    y(n,:) = y(n,:)/y(n,1);
    % residuals = y(n,:) - ySim;
    plot(y(n,1:nPts),'.') 
    % errorbar(y(n,:),repmat(normStdBG,1,24),'.' );  
end

%% Save the figure
ylim([0 1.1])
legend ({'Simulated Signal', 'Sample Pixel 1','Sample Pixel 2','Sample Pixel 3','Sample Pixel 4','Sample Pixel 5','Sample Pixel 6'},'Position',[0.6, 0.7, 0.1,0.1],'location','best')
xlabel 'TE Index'
ylabel 'Normalised Signal'
savefig([savingdir,'/figures/compareSimwithData_Phantom_',phantomName,'__Offset_list_',num2str(offsetListNum),'_compartmentcentercoordslist:',num2str(compartmentCentersList),'.fig'])
matlab2tikz('figurehandle',simCom,'filename',[savingdir,'/DTC_report/',phantomName,'simCom',num2str(offsetListNum),'slice',num2str(sliceNumber)],'height', '\figureheight', 'width', '\figurewidth')
    
%% figure; plot(residuals,'+')

>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end