% This function crops the input image by Crop = [Top, Bottom, Left, Right]
% percent on the each side of the image.  All other dimensions are left
% as before.  If ZeroIm = true then the area outside the crop boundaries
% will be zeroed rather than removed.
%
% function out = CropIm(Im,Crop,ZeroIm);

function out = CropIm(Im,Crop,ZeroIm)
  
  if nargin < 3; ZeroIm = false; end
  
  % Determine image size
  s = size(Im);
  
  % Convert crop margins to position in image
  Crop(2:2:end) = 100 - Crop(2:2:end);
  
  % Convert percentage into pixels (add [1 0 1 0] to prevent reading
  % before the beginning of the array)
  c = round(Crop(1:4)/100 .* [s(1) s(1) s(2) s(2)]) + [ 1 0 1 0];
                                  
  % Crop the array
  if ZeroIm == false
    out = Im( c(1):c(2), c(3):c(4), :, :, : );
  else
    out = zeros(size(Im));
    out(c(1):c(2),c(3):c(4),:,:,:) = Im(c(1):c(2),c(3):c(4),:,:,:);
  end