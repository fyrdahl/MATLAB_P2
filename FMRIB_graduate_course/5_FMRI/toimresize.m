% Resizes an image: alternative to imresize (part of the image processing
% toolbox).
%
function out = toimresize(in,newsize,method)

% Check current size
oldsize = size(in);
xold = 1:oldsize(2);
yold = 1:oldsize(1);

% Define new points for interpolation
xnew = linspace(1,oldsize(2),newsize(2));
ynew = linspace(1,oldsize(1),newsize(1));

% Interpolate
out = zeros([newsize size(in,3)]);
for z = 1:size(in,3)
    out(:,:,z) = interp2(xold,yold,in(:,:,z),xnew,ynew',method);
end