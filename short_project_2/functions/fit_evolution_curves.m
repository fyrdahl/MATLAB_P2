function[varargout] = fit_evolution_curves(phantomName, region, evoType, plotFlag, varargin)
% [varargout] = fit_evolution_curves(phantomName, region, evoType, plotFlag, varargin)
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% 




% Models used to fit the data
T2model = @(a,T2,c, x) a*exp(-x*(1/T2)) + c;
%T1model = @(Mzeq,Mz0,B,T1, x) Mzeq - (Mzeq - Mz0)*B*exp(-x*(1/T1));


switch region
    % fit sample of pixels
    case 'compartments'
        
        switch evoType
            case 'T2'
                % Fit T2 decay
                T2fig = figure;
                for n = 1:6
                    compartmentCenters = varargin{1} ;
                    T2curves(:,n) = squeeze(TEimages(compartmentCenters(n,1,1),compartmentCenters(n,2,1),T2_x));
                    [fittedT2Curve, goodness, output] = fit(T2_x,T2curves(:,n),T2model,'Upper',[5000 4000 500],'Lower',[0 0 -100],'StartPoint',[0.5 200 100]);
                    fittedT2s(n) = fittedT2Curve.T2;
                    varargout{1} = fittedT2s;
                    compartment_a(n) = fittedT2Curve.a;
                    compartment_c(n) = fittedT2Curve.c;
                    
                    switch plotFlag
                        case 'Plot'
                            %Plot and save fits
                            plot(T2_x, T2curves(:,n),'x')
                            hold on
                            filename = '/Users/jallen/Documents/MATLAB/short_project_2/';
                            matlab2tikz('figurehandle',T2fig,'filename',[filename,'/',phantomName,'T2fit',num2str(n)],'height','\figureheight','width','\figurewidth')
                    end
                end
                
                switch plotFlag
                    case 'plot'
                        % Plot and save fits
                        for n = 1:6
                            x = 0:max(T2_x);
                            plot(x, compartment_a(n)*exp(-(x)*(1/fittedT2map(n))) + compartment_c(n),'k-')
                            initialM(n) = compartment_a(n) + compartment_c(n);
                            hold on
                            legend off
                        end
                        legend ({'Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Compartment 5', 'Compartment 6'},'Position',[0.39,0.45,0.25,0.1],'FontSize',8)
                        xlabel 'Echo Time (TE) [ms]'
                        ylabel 'Signal [Arbitrary Units]'
                        matlab2tikz('figurehandle',T2fig,'filename',[filename,'/',phantomName,'T2fitALL'],'height','\figureheight','width','\figurewidth')
                end
                
            case 'T1'
                % Fit T1 decay
                for n = 1:6
                    compartmentCenters = varargin{1} ;
                    fittedT2Curve = 0;
                    goodness = 0;
                    output = 0;
                    T1curves(n,:) = squeeze(TIimages(compartmentCenters(n,1,2),compartmentCenters(n,2,2),T1_x));
                    opts = struct('debug',0,'fiteff',1);
                    [pd(n) r1(n) eff(n) res F(n,:) mz(n,:)] = qmap_t1_fit_ir(T1curves(n,:), 0.001*T1_x,opts);
                    disp(['fit number: ',num2str(n)])
                    fittedT1s(n) = (1/r1(n))*1000
                    varargout{1} = fittedT1s;
                    hold on
                end
                
                switch plotFlag
                    case 'noPlot'
                        %plot T1 fits
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
                        filename = '/Users/jallen/Documents/MATLAB/short_project_2/';
                        legend ({'Compartment 1', 'Compartment 2', 'Compartment 3', 'Compartment 4', 'Compartment 5', 'Compartment 6'},'Location','best')
                        xlabel (['Inversion Time (TI) [(ms)]'])
                        ylabel 'Signal [Arbitrary Units]'
                        matlab2tikz('figurehandle',T1fig,'filename',[filename,'/',phantomName,'T1fitALL'],'height','\figureheight','width','\figurewidth')
                end
        end
        
        %% Fit full phantom, to produce T1 and T2 maps
    case 'fullPhantom'
        
        
        load('mask.mat') %exclude pixels outside of phantom
        mask = reshape(mask,[1,size(mask,1)*size(mask,2)]);
        switch evoType
            case 'T1'
                TIimages = varargin{1};
                T1_x =  varargin{2};
                T1curves = 0;
                T1images = reshape(TIimages,[size(TIimages,1)*size(TIimages,2), size(TIimages,3)]);
                fittedT1map = zeros(1,size(T1images,1));
                fittedBetaMap = zeros(1,size(T1images,1));
                tic
                
                disp(['fitting T1 map...'])
                for r = 1:size(T1images,1)
               % if mask(1,r) > 0
                        %T1 fitting
                        opts = struct('debug',0,'fiteff',1);
                        [pd(r) r1(r) eff] = qmap_t1_fit_ir(T1images(r,T1_x), 0.001*T1_x,opts);
                        fittedBetaMap(1,r) = eff;
                        fittedT1map(1,r) = (1/r1(r))*1000;
             %   end
                end
                
                fittedBetaMap =  reshape(fittedBetaMap, [sqrt(size(fittedBetaMap,2)) , sqrt(size(fittedBetaMap,2)) ]);
                fittedT1map =  reshape(fittedT1map, [sqrt(size(fittedT1map,2)) , sqrt(size(fittedT1map,2)) ]);
                varargout{1} = fittedT1map;
                disp(['Fitting T1 map: complete'])
                toc
        end
        
        
        switch evoType
            
            case 'T2'
                TEimages = varargin{1};
                T2_x = varargin{2};
                T2curves = 0;
                T2images = reshape(TEimages,[size(TEimages,1)*size(TEimages,2), size(TEimages,3)]);
                fittedT2map = zeros(1,size(T2images,1));
                
                tic
                disp('calculating T2 map...')
                for r = 1:size(T2images,1)
                 % if mask(1,r) > 0
                        [fittedT2Curve] = fit(T2_x',T2images(r,T2_x)',T2model,'Upper',[5000 5000 2500],'Lower',[0 0 0],'StartPoint',[1000 100 100]);
                        fittedT2map(1,r) = fittedT2Curve.T2;
                 % end
                end
                toc
                fittedT2map =  reshape(fittedT2map, [sqrt(size(fittedT2map,2)) , sqrt(size(fittedT2map,2)) ]);
                varargout{1} = fittedT2map;
                disp(['Fitting T2 map: complete'])
        end
        
end