function visualise_images(data, range)
% function VISUALISE_IMAGES(DATA,RANGE)
%
% Author: jack.allen@jesus.ox.ac.uk
%
% Description: Function to visualise a series of images, one by one.
%   
% DATA must be 3-dimensional.
%
% RANGE is a 2-element vector, specifying the images to be displayed.
%%

figure;
for i = range(1):range(2)
    if max(max(data(:,:,i))) > 0
        disp(['Index:',num2str(i)])
        colormap gray
        imagesc(data(:,:,i))
        title(['Image Index:',num2str(i)])
        caxis([0 max(max(max(data(:,:,:))))]);
        c = colorbar;
        ylabel(c,'Signal')
        pause
    end
end

end
