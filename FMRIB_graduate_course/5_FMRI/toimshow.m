% Displays an image and returns a handle to the axes
function h = toimshow(img)

h = imagesc(img); axis image; axis off; colormap jet; caxis([0 10]); colorbar;


