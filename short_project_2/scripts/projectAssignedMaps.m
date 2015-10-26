for offsetListNum = 2:8
    %%
    offsetListNum
    clear FA_fig
    clear matchedT1_fig
    clear matchedT2_fig
    clear matchedT1
    clear matchedT2
    clear matchedFAdevInd
    clear M0_mean
    clear scales
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedT1.mat'])
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedT2.mat'])
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'matchedFAdevInd.mat'])
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0_mean.mat'])
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0.mat'])
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'bestMatch.mat'])
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'scales.mat'])
    load([savingdir,'/MAT-files/matches/',phantomName,'list',num2str(offsetListNum),'paramList',num2str(paramList),'M0fit_grad.mat'])
    
    FA_fig = figure
    set(FA_fig,'name',[phantomName,', List',num2str(offsetListNum),', B1 deviation'])
    imagesc(squeeze(matchedFAdevInd(:,:)))
    axis off
    colormap jet
    %   title ([phantomName,', List',num2str(offsetListNum),', B1 deviation'])
    colorbar
    cFA = colorbar;
    cmin = min(dictionaryParams(3,1:sum(dictionaryParams(3,:) ~= 0))) - 0.1*min(dictionaryParams(3,1:(sum(dictionaryParams(3,:) ~= 0))));
    cmax = max(dictionaryParams(3,1:sum(dictionaryParams(3,:) ~= 0)))
    %  cFA.YTick = [cmin , cmax]
    caxis([cmin,cmax])
    ylabel(cFA,'Fraction of B1')
    
    saveas(FA_fig, [savingdir,'/figures/',phantomName,'matchedFAdevIndoffsetlist',num2str(offsetListNum),'.png'])
    matlab2tikz('figurehandle',FA_fig,'filename',[savingdir,'/figures/',phantomName,'slice',num2str(sliceNumber),'FAlist',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    
    pause(2)
    matchedT1_fig = figure
    set(matchedT1_fig,'name',[phantomName,', List',num2str(offsetListNum),', T1'])
    imagesc(matchedT1(:,:))
    axis off
    colormap jet
    %title ([phantomName,', List',num2str(offsetListNum),', T1'])
    cT1 = colorbar;
    switch phantomName
        case 'sphereD170'
            cmin = min(dictionaryParams(1,1:sum(dictionaryParams(1,:) ~= 0))) - 0.1*min(dictionaryParams(1,1:(sum(dictionaryParams(1,:) ~= 0))));
            cmax = max(dictionaryParams(1,1:sum(dictionaryParams(1,:) ~= 0)))
            
        case 'Jack'
            switch ROI
                case 'compartments'
                    cmax = compartmentT1s(2)
                    cmin = 50
                    cmax = 300
                case 'fullPhantom'
                    cmin = 50
                    cmax = 300
            end
    end
    
    cT1.YTick = [cmin : 10 : cmax]
    caxis([cmin,cmax])
    ylabel(cT1,'T1 (ms)')
    
    %saveas(matchedT1_fig, filenameT1 )
    %saveas(matchedT1_fig, [filename,'/',phantomName,'matchedT1offsetlist',num2str(offsetListNum),'-1.png'])
    matlab2tikz('figurehandle',matchedT1_fig,'filename',[savingdir,'/figures/',phantomName,'slice',num2str(sliceNumber),'T1list',num2str(offsetListNum)','ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    
    clear temp
    pause(2)
    matchedT2_fig = figure; imagesc(matchedT2(:,:))
    axis off
    set(matchedT2_fig,'name',[phantomName,', List',num2str(offsetListNum),', T2'])
    colormap jet
    %title ([phantomName,', List',num2str(offsetListNum),', T2'])
    cT2 = colorbar;
    switch phantomName
        case 'sphereD170'
            cmin = min(dictionaryParams(2,1:sum(dictionaryParams(2,:) ~= 0))) - 0.1*min(dictionaryParams(2,1:(sum(dictionaryParams(2,:) ~= 0))));
            cmax = max(dictionaryParams(2,1:sum(dictionaryParams(2,:) ~= 0))) ;
        case 'Jack'
            for i = 1:size(compartmentCenters(:,:,3),1)-1
                temp(i) = matchedT2(squeeze(compartmentCenters(i,1,2)),squeeze(compartmentCenters(i,2,2)))
                % tp(i) = squeeze(compartmentCenters(i,:,3));
            end
            cmin = min(temp)
            cmax = max(temp)
            cmin = 10
            cmax = 110
            
    end
    %   cmax = 120
    % cmax = compartmentT2s(2) + 0.1*compartmentT2s(2)
    cT2.YTick = [cmin : 10 : cmax];
    caxis([cmin,cmax])
    %  caxis([min(dictionaryParams(2,1:sum(dictionaryParams(2,:) ~= 0))) - 0.1*min(dictionaryParams(2,1:(sum(dictionaryParams(2,:) ~= 0)))) ,max(dictionaryParams(2,1:(sum(dictionaryParams(2,:) ~= 0)))) ])
    ylabel(cT2,'T2 (ms)')
    %saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum)])
    saveas(matchedT2_fig, [savingdir,'/figures/',phantomName,'matchedT2offsetlist',num2str(offsetListNum),'.png'])
    matlab2tikz('figurehandle',matchedT2_fig,'filename',[savingdir,'/figures/',phantomName,'slice',num2str(sliceNumber),'T2list',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    pause(2)
    
    image = log10(abs(M0_mean(:,:)));
    M0_mean_fig = figure; imagesc(image(:,:))
    axis off
    set(M0_mean_fig,'name',[phantomName,', List',num2str(offsetListNum),', M0_mean'])
    colormap jet
    cM0_mean = colorbar;
    %    for i = 1:size(compartmentCenters(:,:,3),1)
    %        temp(i) = image(squeeze(compartmentCenters(i,1,3)),squeeze(compartmentCenters(i,2,3)));
    %         tp(i) = squeeze(compartmentCenters(i,:,3));
    %   end
    %   cmin = min(temp)
    %   cmax = max(temp)
    %cmax = 120
    %caxis([0,4500])
    ylabel(cM0_mean,'M0 (log_{10}(Mean of Absolute Scaling Factors))')
    saveas(M0_mean_fig, [savingdir,'/figures/',phantomName,'matchedM0_mean',num2str(offsetListNum),'.png'])
    matlab2tikz('figurehandle',M0_mean_fig,'filename',[savingdir,'/figures/',phantomName,'slice',num2str(sliceNumber),'M0mean',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    
    %     pause(2)
    %     M0fig = figure; imagesc(log10(M0(:,:)))
    %     axis off
    %     set(M0fig,'name',[phantomName,', List',num2str(offsetListNum),', M0'])
    %     colormap jet
    %     cT2 = colorbar;
    %     ylabel(cT2,'log_{10}(M0)')
    %
    %     %saveas(matchedT2_fig, [workingdir,'/figures/',phantomName,'matchedT2_offsetList',num2str(offsetListNum)])
    %     saveas(M0fig, [filename,'/',phantomName,'matchedM0',num2str(offsetListNum),'.png'])
    %     matlab2tikz('figurehandle',M0fig,'filename',[filename,'/',phantomName,'M0',num2str(offsetListNum)],'height', '\figureheight', 'width', '\figurewidth')
    %
    
    
    data = reshape(data,[64*64, 24]);
    bestMatch(1001,:);
    M0fit_grad_fig = figure; imagesc(M0fit_grad(:,:));
    axis off
    set(M0fit_grad_fig,'name',[phantomName,', List',num2str(offsetListNum),', M0fit_grad'])
    colormap jet
    cM0fit_grad = colorbar;
    ylabel(cM0fit_grad,'M_{0}R [a.u.]')
    saveas(M0fit_grad_fig, [savingdir,'/figures/',phantomName,'M0fit_grad',num2str(offsetListNum),'.png'])
    matlab2tikz('figurehandle',M0fit_grad_fig,'filename',[savingdir,'/figures/',phantomName,'slice',num2str(sliceNumber),'M0fit_grad',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    
end