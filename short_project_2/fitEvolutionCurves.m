function[compartmentT1s, compartmentT2s, fittedCurve, goodness, output] = fitEvolutionCurves(TEimages, TIimages, T2_x, T1_x, region, varargin)

% Models to fit to the data
% T2model = @(a,T2,c, x) a*exp(-x*(1/T2)) + c
T2model = @(a,T2,c, x) a + -x*(1/T2) + c
T1model = @(Mzeq,Mz0,B,T1, x) Mzeq - (Mzeq - Mz0)*B*exp(-x*(1/T1))

switch region
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'compartments'
        
        %% Fit T2 decay
        for n = 1:6
            compartmentCenters = varargin{1} ;
            %subtract noise
            curve2 = log( squeeze(TEimages(compartmentCenters(n,1),compartmentCenters(n,2),T2_x)))
            %curve = squeeze(mean(mean(images(23:47,22:43,x))));
            
            [fittedCurve, goodness, output] = fit(T2_x,curve2,T2model,'Upper',[5000 4000 500],'Lower',[0 0 -100],'StartPoint',[0.5 200 100])
            
            compartmentT2s(n) = fittedCurve.T2
            
            figure; plot(T2_x, curve2,'.')
            hold on
            plot(fittedCurve)
        end
        
        %% Fit T1 decay
        for n = 1:6
            compartmentCenters = varargin{1} ;
            curve1 = squeeze(TIimages(compartmentCenters(n,1),compartmentCenters(n,2),T1_x));
            %curve = [compartmentTIimages(n,:)]'
            %                 [fittedCurve, goodness, output] = fit(x,curve,T1model,'Upper',[2000 2000 1 5000],'Lower',[0 -5000 0.9 0],'StartPoint',[1000 100 0.5 100]);
            %                 compartmentT1s(n) = fittedCurve.T1
            
            fittedCurve = 0;
            goodness = 0;
            output = 0;
            [pd r1 eff res] = shurleyT1fit(curve1', 'compartments', 1, T1_x);
            disp(['fit number: ',num2str(n)])
            compartmentT1s(n) = (1/r1)*1000
            
            %                 figure; plot(x, curve,'.')
            %                 hold on
            %                 plot(fittedCurve)
        end
        
        
        
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