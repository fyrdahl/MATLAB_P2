function [] = show_pair(im1,im2, fov_square)

% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

% show_pair(im1, im2)
%
% display images im1 and im2 next to each other (im1 on left)
% fov_square is an optional boolean specifying whether to fix the FOV of 
% the second subplot to be square

if nargin < 3
    fov_square = 0;
end


     % display data in both domains
     figure;
     subplot(121); show_im(abs(im1));
     subplot(122); show_im(abs(im2));
     if fov_square
         axis square;
     end


function [l] = show_im(im)

% show_im(im)
%
% displays an image

     % display contrast
     m = mean(im(:)); s = std(im(:));
     low  = max(m-2*s, min(im(:)));
     high = min(m+2*s, max(im(:)));

     % grayscale colormap 
     a=[1 1 1]/256; b=[1:256]; c=(a'*b)';
     colormap('gray');

     % display image
     if (isreal(im))
       imagesc(rot90(im),[low,high]);
     else
       imagesc(rot90(abs(im)),[low,high]);
     end;
     set(gca,'YDir','normal');
     axis('equal'); axis([1 size(im,1) 1 size(im,2)]);

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')
