function plot_sim_comparison_Hargreaves(data, T1, T2, offsets,nRuns,df)
% plot_sim_comparison(Mxy,data,nPts,phantomName,savingdir,offsetListNum, sliceNumber)
%


%% simulate the signal using Bernstein mathematical description
TRmin = 130;
TEmin = 32;
%[Msig, M] = sim_SE_Hargreaves(T1, T2, minTR, minTE, offsets,nRuns, df, plotFlag)
Msig = sim_SE_Hargreaves(T1, T2, TRmin, TEmin, offsets,nRuns, df,'noPlot');
% normalise to the first image index
simSig = Msig(:)/Msig(1);


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
    y = data(n,m,:);
    y = y(:)/y(1);
    plot(y,'-o')
    hold on
   end
end

    % errorbar(y(n,:),repmat(normStdBG,1,24),'.' );
%end
plot(simSig,'b-*','LineWidth',3,'MarkerSize',15)
hold on

ylim([0 max(max(y))])
%legend ({'Simulated Signal', 'Sample Pixel 1','Sample Pixel 2','Sample Pixel 3','Sample Pixel 4','Sample Pixel 5','Sample Pixel 6'},'location','best')
xlabel 'TE Index'
ylabel 'Normalised Signal'

set(gca, 'FontSize',18)

%%
%savefig([savingdir,'/figures/compareSimwithData_Phantom_',phantomName,'__Offset_list_',num2str(offsetListNum),'_compartmentcentercoordslist:',num2str(compartmentCentersList),'.fig'])
%matlab2tikz('figurehandle',simCom,'filename',[savingdir,'/DTC_report/',phantomName,'simCom',num2str(offsetListNum),'slice',num2str(sliceNumber)],'height', '\figureheight', 'width', '\figurewidth')
% figure; plot(residuals,'+')
end