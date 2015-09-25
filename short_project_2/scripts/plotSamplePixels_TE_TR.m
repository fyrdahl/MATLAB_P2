figure;
imagesc(TEimages(:,:,TE(2)))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1 :plotNumCompartments
    plot(compartmentCenters(i,2,1),compartmentCenters(i,1,1),'*')
    text(compartmentCenters(i,2,1),compartmentCenters(i,1,1), compartmentLabels(i) )
end

TIfig = figure;
imagesc(TIimages(:,:,TI(2)))
hold on
compartmentLabels = ['1', '2','3','4','5','6'];
for i = 1:plotNumCompartments
    plot(compartmentCenters(i,2,2),compartmentCenters(i,1,2),'*')
    t = text(compartmentCenters(i,2,2),compartmentCenters(i,1,2), compartmentLabels(i))
    switch phantomName
        case 'Jack'
      set(t,'color','white')
        case 'sphereD170'            
    set(t,'color','k')
    end
end
axis off
colormap gray
set(TIfig,'name','TI sample pixels')
figFilepath = '/Users/jallen/Documents/MATLAB/short_project_2/DTC_report'
matlab2tikz('figurehandle',TIfig,'filename',[figFilepath,'/',phantomName,'TILabels'],'height','\figureheight','width','\figurewidth')