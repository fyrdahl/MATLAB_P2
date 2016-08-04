function make_TC_video(data,videoName,varargin)

% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
%
% Make video from TC images.


v = VideoWriter(videoName);
if nargin ==2
v.FrameRate = 2;
else
    v.FrameRate = varargin{1};
end
open(v)

for k = 1:size(data,3) 

    title(['Frame ',num2str(k)])
   imagesc(data(:,:,k))
   colormap gray
   colorbar
   max(max(data(:,:,1)))
   caxis([ 0 max(max(data(:,:,1)))])
   pause(0.5)
   drawnow
   frame = getframe;
   writeVideo(v,frame);
end

close(v);

end
