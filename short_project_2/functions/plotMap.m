function plotMap(phantomName,mapType, offsetListNum, varargin)

map = load([phantomName,'list',num2str(offsetListNum),'matched',mapType,'.mat']);

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


switch mapType 
    case 'FAdevInd'
     caxis([ 0.5 1.5])
end
    
set(fig,'name',[phantomName,', List',num2str(offsetListNum),' ',mapType],'numbertitle','off')

colorbar
axis off

%varagin{2} is the working directory
saveas(fig, [varargin{1},'/figures/',phantomName,'_list',num2str(offsetListNum),'_',mapType,'_map.png'])

end