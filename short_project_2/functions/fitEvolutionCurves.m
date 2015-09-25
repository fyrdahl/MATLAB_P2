function[compartmentT1s, compartment_T2s, T2curves, T1curves, fittedCurve, goodness, output,F] = fitEvolutionCurves(phantomName,TEimages, TIimages, T2_x, T1_x, region, varargin)

% Models to fit to the data
T2model = @(a,T2,c, x) a*exp(-x*(1/T2)) + c;
%T2model = @(a,T2,c, x) a + -x*(1/T2);
%T1model = @(Mzeq,Mz0,B,T1, x) Mzeq - (Mzeq - Mz0)*B*exp(-x*(1/T1));

switch region
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'compartments'
        
        %% Fit T2 decay
        T2fig = figure;
        for n = 1:6
            compartmentCenters = varargin{1} ;
            
            T2curves(:,n) = squeeze(TEimages(compartmentCenters(n,1,1),compartmentCenters(n,2,1),T2_x));
            
            [fittedCurve, goodness, output] = fit(T2_x,T2curves(:,n),T2model,'Upper',[5000 4000 500],'Lower',[0 0 -100],'StartPoint',[0.5 200 100]);
            
            compartment_T2s(n) = fittedCurve.T2;
            compartment_a(n) = fittedCurve.a;
            compartment_c(n) = fittedCurve.c;
            
            plot(T2_x, T2curves(:,n),'.')
            
            hold on
            
            filename = '/Users/jallen/Documents/MATLAB/short_project_2/DTC_report'
            matlab2tikz('figurehandle',T2fig,'filename',[filename,'/',phantomName,'T2fit',num2str(n)],'height','\figureheight','width','\figurewidth')
        end
        
        for n = 1:6
            plot(T2_x, compartment_a(n)*exp(-(T2_x)*(1/compartment_T2s(n))) + compartment_c(n),'k-+')
            hold on
            legend off
        end
        legend ({'Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Compartment 5', 'Compartment 6'},'Position',[0.39,0.45,0.25,0.1],'FontSize',8)
        
        xlabel 'Echo Time (TE) [ms]'
        ylabel 'Signal [Arbitrary Units]'
        
        matlab2tikz('figurehandle',T2fig,'filename',[filename,'/',phantomName,'T2fitALL'],'height','\figureheight','width','\figurewidth')
        
        
        clear T2fig
        %% Fit T1 decay
        figure
        for n = 1:6
            
            compartmentCenters = varargin{1} ;
            
            fittedCurve = 0;
            goodness = 0;
            output = 0;
            
            T1curves(n,:) = squeeze(TIimages(compartmentCenters(n,1,2),compartmentCenters(n,2,2),T1_x));
            opts = struct('debug',1,'fiteff',1);
            
            [pd(n) r1(n) eff res F(n,:) mz(n,:)] = qmap_t1_fit_ir(T1curves(n,:), 0.001*T1_x,opts);
            
            disp(['fit number: ',num2str(n)])
            compartmentT1s(n) = (1/r1(n))*1000
            hold on
            
        end
        
        T1fig = figure;
        for n = 1:6
            plot(log10(T1_x), T1curves(n,:) ,'.')
            hold on
        end
        
        for n = 1:6
            r1(n)
            pd(n)
            T1 = (1/r1(n))*1000;
            x = 1:numel(F(1,:));
            plot(log10(T1_x), F(n,x) ,'k+' )
            plot(log10(T1_x), mz(n,x) ,'k-*' )
            hold on
        end
        legend off
        
        
        %%
        filename = '/Users/jallen/Documents/MATLAB/short_project_2/DTC_report'
        
        
        legend ({'Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Compartment 5', 'Compartment 6'},'Location','best')
        
        xlabel (['Inversion Time (TI) [log_{10}','(ms)]'])
        ylabel 'Signal [Arbitrary Units]'
       ylim([ -2000 3500])
        matlab2tikz('figurehandle',T1fig,'filename',[filename,'/',phantomName,'T1fitALL'],'height','\figureheight','width','\figurewidth')
        
        %
        %% Fit full phantom, to produce T1 and T2 maps
        
    case 'fullPhantom'
        load('mask.mat')
        
        for r = 1:64
            disp('calculating T2 map...')
            disp(['row ',num2str(r),' of 64'])
            for c = 1:64
                if mask(r,c) >0
                    
                    curve = squeeze(images(r,c,T1_x));
                    
                    [fittedCurve, goodness, output] = fit(T1_x',curve,T2model,'Upper',[5000 5000 2500],'Lower',[0 0 0],'StartPoint',[1000 100 100]);
                    
                    T2fits(r,c,1) = fittedCurve.a;
                    T2fits(r,c,2) = fittedCurve.T2;
                    T2fits(r,c,3) = fittedCurve.c;
                    T2fitResiduals(r,c,:) = output.residuals;
                    T2map(r,c) = fittedCurve.T2;
                    
                    curve = squeeze(images(r,c,T2_x));
                    %                     [fittedCurve, goodness, output] = fit(x,curve,T1model,'Upper',[2000 2000 10 5000],'Lower',[-5000 -5000 0 0],'StartPoint',[1000 100 0.5 100])
                    %                     T1map(r,c) = fittedCurve.T1;
                    [pd r1 eff res] = shurleyT1fit(curve', 'compartments', 1, TI);
                    T1map(r,c) = 1./r1;
                    
                end
                
            end
        end
        disp('finished calculating T2 map')
        toc
        
        if relaxationType == 'T2' && figureFlag == 'showMap'
            figure; imagesc(T2map(:,:))
            caxis([0 150])
            
        end
        
        if relaxationType == 'T1' && figureFlag == 'showMap'
            figure; imagesc(T1map(:,:))
            caxis([0 350])
        end
        
end

end