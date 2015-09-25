function [pd r1 eff res] = shurleyT1fit(images, plotFlag,TI)
%
% function to use Sam Hurley's qmap code

        for i = 2:numel(TI)
            data = reshape(imagesn[size(images;
            images(:,i) = d(:);
        end
        
        opts = struct('fiteff',1)
        
        if plotFlag == 'plot'
            clear opts
            opts = struct('debug',1,'fiteff',1)
        end
        
        [pd r1 eff res] = qmap_t1_fit_ir(images, 0.001*TI(2:size(images,2)),opts)
        

end

