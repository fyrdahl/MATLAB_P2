<<<<<<< HEAD
function visualise_images(data, range)
% function VISUALISE_IMAGES(DATA,RANGE)
=======
function visualise_Images(data, range)
% function visualise_Images(DATA)
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
%
% Author: jack.allen@jesus.ox.ac.uk
%
% Description: Function to visualise a series of images, one by one.
<<<<<<< HEAD
%   
% DATA must be 3-dimensional.
%
% RANGE is a 2-element vector, specifying the images to be displayed.
=======
% DATA must be 3-dimensional
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
%%

figure;
for i = range(1):range(2)
    if max(max(data(:,:,i))) > 0
<<<<<<< HEAD
        disp(['Index:',num2str(i)])
=======
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
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
