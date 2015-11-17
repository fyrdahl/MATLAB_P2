function visualiseImages(data)
%% Jack Allen
% jack.allen@jesus.ox.ac.uk
% Function to visualise a series of images, one by one.
%
% DATA must be 3-dimensional
%%

figure;
for i = 1:size(data,3)
  
    if max(max(data(:,:,i))) > 0
 colormap gray
 title ([num2str(i)])
imagesc(data(:,:,i))
caxis([0 max(max(max(data(:,:,:))))]);
colorbar
pause
    end
end

end
