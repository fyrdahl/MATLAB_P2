% Creates distortion and dropout in the input image based on the fieldmap
% and other parameters, for EPI readouts.  Based on code by James Meakin,
% 2011.
%
% Tom Okell, 2011: tokell@fmrib.ox.ac.uk


% Wenchuan Wu, 2014 modified for dwi tpi

function out = Distortion_dwiEPI(inallo, FMap, UpSampFactor, TE, BWPerPix, MtxSize, PhaseEncDir, blip_dir)


no_dir = size(inallo,3);

if strcmp(PhaseEncDir,'LR')
    
    inall = zeros(size(inallo,2),size(inallo,1),size(inallo,3));
    
    for ii = 1 : no_dir
   inall(:,:,ii)  = rot90(inallo(:,:,ii));
   
   
   MtxSize = size(inall);
    end
     FMap = rot90(FMap); 
     inallo = inall;
else
   
    inall = inallo;
    
end
    


Phase_rads = TE/1000 * repmat(FMap,[1,1,no_dir]) * 2*pi;
inall = inall .* exp(1i*Phase_rads);

if UpSampFactor ~=1;
  Sz = size(inall);
  if length(Sz) == 2
     Sz = [Sz 1]; 
  end
  
  
  [x,y] = meshgrid(1:Sz(2),1:Sz(1));
  dxnew = 1/UpSampFactor;
  [xnew, ynew] = meshgrid((0.5+dxnew/2):dxnew:(Sz(2)+0.5-dxnew/2), ...
                          (0.5+dxnew/2):dxnew:(Sz(1)+0.5-dxnew/2));
  
  inall = zeros(size(xnew,1),size(xnew,2),Sz(3));

  for dd = 1 : Sz(3)

  inall(:,:,dd)   = interp2(x,y,inallo(:,:,dd),  xnew,ynew,'*spline'); % faster hopefully

  end
                      
  FMap = interp2(x,y,FMap,xnew,ynew,'*spline');
     inall = inall/UpSampFactor^2;

end

FMap(find(isnan(FMap))) = 0;

% Record the new size
Sz = size(inall);

%% end






% Get the positions of the voxel
[x,y] = meshgrid(1:Sz(2),1:Sz(1));

% Initalise output
% out = zeros(round(size(in)/max(Sz) * MtxSize));

SzOut = [MtxSize(1) MtxSize(2) no_dir];
out = zeros(SzOut);

% Calculate the effective bandwidth in the phase encode direction
BWPerPixPE = BWPerPix / SzOut(2);

% Using the field map, find the new position of each voxel
if strcmp(blip_dir,'Down')
    Factor = -1; 
else
    Factor = 1;
end
ynew = y;
xnew = x + Factor * FMap / BWPerPixPE * Sz(2) / SzOut(2);

% Account for aliasing
xnew = mod(xnew,Sz(2));

% Sum/mean in the readout direction (same index for each column)
dx = Sz(2)/SzOut(2);           dy = Sz(1)/SzOut(1);
x_bounds = 0.5:dx:(Sz(2)+0.5); y_bounds = 0.5:dy:(Sz(1)+0.5);


% toc
% tic

tmp = zeros(SzOut(1),Sz(2),no_dir); tmp_x = zeros(SzOut(1),Sz(2));

% for dd = 1:SzOut(3)  % avoid using for loop
for jj = 1:SzOut(1)
    Idx = ( ynew(:,1) >= y_bounds(jj) ) & ( ynew(:,1) < y_bounds(jj+1) );
    

      tmp(jj,:,:) = sum(inall(Idx,:,:),1);
    tmp_x(jj,:) = mean(xnew(Idx,:),1);
end
% end
inall = tmp; xnew = tmp_x;



for ii = 1:SzOut(2) % Loop over x
  for jj = 1:SzOut(1) % Loop over y
    Idx = ( xnew(jj,:) >= x_bounds(ii) ) & ( xnew(jj,:) < x_bounds(ii+1) );
    
      
      out(jj,ii,:) = abs(sum(inall(jj,Idx,:),2));  % avoid using for loop
      

  end
end
% end

% Unrotate the output
if strcmp(PhaseEncDir,'LR')
    
    out2 = zeros(size(out,2),size(out,1),size(out,3));
    for dd = 1 : no_dir
   out2(:,:,dd)   = rot90(out(:,:,dd),-1);
    end
    out = out2;
end
% toc
% Taking the mean of no numbers leads to NaN so replace with zeros
out(isnan(out)) = 0;