function plotMap(phantomName,mapType, offsetListNum, varargin)

map = load(['list',num2str(offsetListNum),'matched',mapType,'.mat']);

img = getfield(map,['matched',mapType]);
fig = figure;
colormap jet
imagesc(img)
           
    
if nargin >= 4
% caxis([varargin{2}])
% else
compartmentCenters = varargin{2}
    for r = 1:size(compartmentCenters,1)-1
    
        samples(r) = img(compartmentCenters(r,1,3), compartmentCenters(r,2,3));

end
    caxis([min(samples)-(max(samples)*0.1), max(samples) + (max(samples)*0.1)])
end
title (['offset list',num2str(offsetListNum),', ',mapType,' map'])

%varagin{2} is the working directory
saveas(fig, [varargin{1},'/figures/',phantomName,'_list',num2str(offsetListNum),'_',mapType,'_map.png'])
end