function [pd r1 eff res] = shurleyT1fit(images, region, plotFlag,TI)

switch region
    case 'fullImage'
        for i = 2:numel(TI)
            d = images(:,:,TI(i));
            data(:,i) = d(:);
        end
        
        if plotFlag == 1
            opts = struct('debug',1)
        end
        
        [pd r1 eff res] = qmap_t1_fit_ir(data, 0.001*TI(2:size(data,2)),opts)
        



case 'compartments'
    for n = 1:size(images,1)
        data(n,:) = images(n,:)
    end
    
    if plotFlag == 1
        opts = struct('debug',1)
    end
    
    [pd r1 eff res] = qmap_t1_fit_ir(data, 0.001*TI(2:size(data,2)),opts)
    
end

