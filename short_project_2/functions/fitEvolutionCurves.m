function[fittedT1s, fittedT2s, T2curves, T1curves, fittedCurve, goodness, output,F] = fit_evolution_curves(phantomName,TEimages, TIimages, T2_x, T1_x, region, varargin)

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
            
            fittedT2s(n) = fittedCurve.T2;
            compartment_a(n) = fittedCurve.a;
           compartment_c(n) = fittedCurve.c;
            
            plot(T2_x, T2curves(:,n),'x')
            
            hold on
            
            filename = '/Users/jallen/Documents/MATLAB/short_project_2/';
            matlab2tikz('figurehandle',T2fig,'filename',[filename,'/',phantomName,'T2fit',num2str(n)],'height','\figureheight','width','\figurewidth')
        end
        
        for n = 1:6
            x = 0:max(T2_x);
           plot(x, compartment_a(n)*exp(-(x)*(1/fittedT2s(n))) + compartment_c(n),'k-')
           initialM(n) = compartment_a(n) + compartment_c(n);
            hold on
            legend off
        end
        legend ({'Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Compartment 5', 'Compartment 6'},'Position',[0.39,0.45,0.25,0.1],'FontSize',8)
        
        xlabel 'Echo Time (TE) [ms]'
        ylabel 'Signal [Arbitrary Units]'
        
        matlab2tikz('figurehandle',T2fig,'filename',[filename,'/',phantomName,'T2fitALL'],'height','\figureheight','width','\figurewidth')

        %% Fit T1 decay
        for n = 1:6
            
            compartmentCenters = varargin{1} ;
            
            fittedCurve = 0;
            goodness = 0;
            output = 0;
            
            T1curves(n,:) = squeeze(TIimages(compartmentCenters(n,1,2),compartmentCenters(n,2,2),T1_x));
            opts = struct('debug',0,'fiteff',1);
            
            [pd(n) r1(n) eff(n) res F(n,:) mz(n,:)] = qmap_t1_fit_ir(T1curves(n,:), 0.001*T1_x,opts);
            
            disp(['fit number: ',num2str(n)])
            fittedT1s(n) = (1/r1(n))*1000
            hold on
            
        end
        
        figure
        for n = 1:6
            plot((T1_x), T1curves(n,:) ,'x')
            hold on
        end
        x = 0:max(T1_x)+0.2*max(T1_x);
        for n = 1:6
            T1(n) = (1/r1(n))*1000;
            y(n,:) = pd(n) - (2*pd(n))*eff(n)*exp(-x*(1/T1(n)));
            plot(x,y(n,:));
            minY(n) = y(n,1);
            maxY(n) = y(n,end);
            hold on
        end
        ylim([min(minY)+0.2*min(minY),max(maxY) + 0.2*max(maxY)])
        xlim([0 , max(x)])
        
        
        %%
        filename = '/Users/jallen/Documents/MATLAB/short_project_2/'
        
        
        legend ({'Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Compartment 5', 'Compartment 6'},'Location','best')
        
        xlabel (['Inversion Time (TI) [(ms)]'])
        ylabel 'Signal [Arbitrary Units]'
        matlab2tikz('figurehandle',T1fig,'filename',[filename,'/',phantomName,'T1fitALL'],'height','\figureheight','width','\figurewidth')
        
        %
        %% Fit full phantom, to produce T1 and T2 maps
        
    case 'fullPhantom'
    
        load('mask.mat')
        T2curves = 0;
        T1curves = 0;
        
        T2images = reshape(TEimages,[size(TEimages,1)*size(TEimages,2), size(TEimages,3)]);
        T1images = reshape(TIimages,[size(TIimages,1)*size(TIimages,2), size(TIimages,3)]);
        mask = reshape(mask,[1,size(mask,1)*size(mask,2)]);
        fittedT1s = zeros(1,size(T1images,1))
        fittedT2s = zeros(1,size(T2images,1))
        fittedBeta = zeros(1,size(T1images,1))
        tic
        for r = 1:size(T2images,1)
            disp('calculating T2 map...')
            disp(['row ',num2str(r),' of ',num2str(size(mask,2)) ])
            
            if mask(1,r) > 0
                
                %T2
                [fittedCurve, goodness, output] = fit(T2_x,T2images(r,T2_x)',T2model,'Upper',[5000 5000 2500],'Lower',[0 0 0],'StartPoint',[1000 100 100]);
                
               % T2fits(r,1) = fittedCurve.a;
              % T2fits(r,2) = fittedCurve.T2;
              %  T2fits(r,3) = fittedCurve.c;
              %  T2fitResiduals(r,:) = output.residuals;
                fittedT2s(1,r) = fittedCurve.T2;
                
                %T1
                opts = struct('debug',0,'fiteff',1);
                [pd(r) r1(r) eff res F(r,:) mz(r,:)] = qmap_t1_fit_ir(T1images(r,T1_x), 0.001*T1_x,opts); 
                 fittedBeta(1,r) = eff;
                 fittedT1s(1,r) = (1/r1(r))*1000;
                

            end
            
        end
        toc
      B1map =  reshape(fittedBeta, [sqrt(size(fittedBeta,2)) , sqrt(size(fittedBeta,2)) ]);
      T2map =  reshape(fittedT2s, [sqrt(size(fittedT2s,2)) , sqrt(size(fittedT2s,2)) ]);
       T1map =  reshape(fittedT1s, [sqrt(size(fittedT1s,2)) , sqrt(size(fittedT1s,2)) ]);
        %  if relaxationType == 'T2' && figureFlag == 'showMap'
        T2fig = figure; imagesc(T2map(:,:))
        caxis([180 300])
     %   ylabel(T2fig,'Gold Standard T2')
        
        %  end
        
        %  if relaxationType == 'T1' && figureFlag == 'showMap'
       T1fig = figure; imagesc(T1map(:,:))
       caxis([180 300])
 %      ylabel(T1fig,'Gold Standard T1')
        %  end
        B1fig = figure; imagesc(B1map(:,:))
 %      ylabel(B1map,'Gold Standard T1')
        
          fittedCurve = 0;
        goodness = 0;
        output = 0;
        F = 0;
end

end