if ROI == 'fullPhantom'
    goldStdT1maps = reshape(compartmentT1s,[sqrt(size(compartmentT1s,2)), sqrt(size(compartmentT1s,2))]);
    goldStdT2maps = reshape(compartmentT2s,[sqrt(size(compartmentT2s,2)), sqrt(size(compartmentT2s,2))]);
    save([savingdir,'/MAT-files/images/',phantomName,'goldStd','T1.mat'],'goldStdT1maps')
    save([savingdir,'/MAT-files/images/',phantomName,'goldStd','T2.mat'],'goldStdT2maps')
    
    clear goldStdT1maps
    clear goldStdT2maps
    load([savingdir,'/MAT-files/images/',phantomName,'goldStd','T1.mat'])
    load([savingdir,'/MAT-files/images/',phantomName,'goldStd','T2.mat'])
    
    filename = '/Users/jallen/Documents/MATLAB/short_project_2/DTC_report/';
    
    switch phantomName
        case 'Jack'
            goldStdT1mapsFig = figure; imagesc(goldStdT1maps)
            axis off
            colormap jet
            cgoldStdT1mapsFig = colorbar;
            ylabel(cgoldStdT1mapsFig,'T1 [ms]')
            cmin = 50;
            cmax = 300;
            caxis([cmin cmax])
            cgoldStdT1mapsFig.YTick = [cmin : 50 : cmax];
            matlab2tikz('figurehandle',goldStdT1mapsFig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'goldStdT1maps',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
            
            goldStdT2mapsFig = figure; imagesc(goldStdT2maps)
            axis off
            colormap jet
            cgoldStdT2mapsFig = colorbar;
            ylabel(cgoldStdT2mapsFig,'T2 [ms]')
            cmin = 10;
            cmax = 110;
            caxis([cmin cmax])
            cgoldStdT2mapsFig.YTick = [cmin : 20 : cmax];
            matlab2tikz('figurehandle',goldStdT2mapsFig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'goldStdT2maps',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
        case 'sphereD170'
            goldStdT1mapsFig = figure; imagesc(goldStdT1maps) 
            axis off
            colormap jet
            cgoldStdT1mapsFig = colorbar;
            ylabel(cgoldStdT1mapsFig,'T1 [ms]')
            cmin = 180;
            cmax = 300;
            caxis([cmin cmax])
            cgoldStdT1mapsFig.YTick = [cmin : 20 : cmax];
            matlab2tikz('figurehandle',goldStdT1mapsFig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'goldStdT1maps',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
            
            goldStdT2mapsFig = figure; imagesc(goldStdT2maps)
            axis off
            colormap jet
            cgoldStdT2mapsFig = colorbar;
            ylabel(cgoldStdT2mapsFig,'T2 [ms]')
            caxis([cmin cmax])
            cgoldStdT2mapsFig.YTick = [cmin : 20 : cmax];
            matlab2tikz('figurehandle',goldStdT2mapsFig,'filename',[filename,'/',phantomName,'slice',num2str(sliceNumber),'goldStdT2maps',num2str(offsetListNum),'ParamList',num2str(paramList)],'height', '\figureheight', 'width', '\figurewidth')
    end
    
end