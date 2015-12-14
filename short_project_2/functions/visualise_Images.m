function visualise_Images(data, range)
% function visualise_Images(DATA)
%
% Author: jack.allen@jesus.ox.ac.uk
%
% Description: Function to visualise a series of images, one by one.
% DATA must be 3-dimensional
%%

figure;
for i = range(1):range(2)
    if max(max(data(:,:,i))) > 0
        colormap gray
        imagesc(data(:,:,i))
        title ([num2str(i)])
        caxis([0 max(max(max(data(:,:,:))))]);
        c = colorbar;
        ylabel(c,'Signal')
        pause
    end
end






end
