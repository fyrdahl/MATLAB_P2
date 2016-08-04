function [mask] = make_roi(image,title_text,max_intensity)
% NNG, 25/10/2013
% replacement for mask_roi for EP practical (where roipoly is not avaliable
% due to lack of image processing toolbox). 
% MChiew - Completed function with patch and inpolygon
% T.Okell - allowed max intensity to be modified

if nargin < 2; title_text = 'Draw polygon to generate mask'; end
if nargin < 3; max_intensity = max(image(:)); end

display('Draw a polygon to make a mask by clicking on the image (first line will appear after second click)');  
display('To finish click the first point again or simply double click'); 

scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]); 
imagesc(image,[0 max_intensity]); colormap('gray'); title(title_text); 
axis image;
hold on;
[x,y] = ginput(2);
plot([x(1),x(2)],[y(1),y(2)],'Color','r','LineWidth',1);

i = 2; 
while ((x(i)~= x(i-1))||(y(i)~= y(i-1))) % finish when start point is clicked again
    i = i + 1; 
    [x(i), y(i)] = ginput(1); 
    plot([x(i-1),x(i)],[y(i-1),y(i)],'Color','r','LineWidth',1);
end  

patch(x,y,'r','FaceAlpha',0.7);

pause(1);
close;

xidx = repmat([1:size(image,1)]',[1,size(image,2)]);
yidx = repmat([1:size(image,2)],[size(image,1),1]);
mask = inpolygon(yidx, xidx, x, y);

%inpolygon
  
%close;
