% Creates distortion and dropout in the input image based on the fieldmap
% and other parameters, for EPI readouts.  Based on code by James Meakin,
% 2011.
%
% Tom Okell, 2011: tokell@fmrib.ox.ac.uk

function out = DistortionAndDropout(in, FMap, UpSampFactor, TE, BWPerPix, MtxSize, SlcThk, PhaseEncDir, blip_down, IsSignal, B0, DropOutMap)


% Through slice dephasing
if IsSignal
  in = in .* exp(-DropOutMap*(B0/7)^2*TE/100*SlcThk/10*50);
  %DispIm([in DropOutMap]);
end

% Rotate the image if necessary so that readout is along rows
if strcmp(PhaseEncDir,'AP')
   in   = rot90(in);
   FMap = rot90(FMap);
end

% Upsample the input image and field map
if UpSampFactor ~=1;
  Sz = size(in);
  [x,y] = meshgrid(1:Sz(2),1:Sz(1));
  dxnew = 1/UpSampFactor;
  [xnew, ynew] = meshgrid((0.5+dxnew/2):dxnew:(Sz(2)+0.5-dxnew/2), ...
                          (0.5+dxnew/2):dxnew:(Sz(1)+0.5-dxnew/2)  );
  
  in   = interp2(x,y,in,  xnew,ynew);
  FMap = interp2(x,y,FMap,xnew,ynew);
  
  % If this is signal, reduce in proportion to the new voxel size
  if IsSignal
      in = in/UpSampFactor^2;
  end
  
end

% Record the new size
Sz = size(in);

% Calculate the phase accrual at each small voxel and apply to the input
% image (only for signals, not parameter maps)
if IsSignal
  Phase_rads = TE/1000 * FMap * 2*pi;
  in = in .* exp(1i*Phase_rads);
  % in = in .* exp(1i*Phase_rads*SlcThk/10); % Fudge for through plane dephasing
  % in = in .* exp(1i*Phase_rads.^2/(6*2*pi)*SlcThk/5.*rand(size(Phase_rads))*B0/3); % New fudge for through plane dephasing  
end

% Get the positions of the voxel
[x,y] = meshgrid(1:Sz(2),1:Sz(1));

% Initalise output
out = zeros(round(size(in)/max(Sz) * MtxSize));
SzOut = size(out);

% Calculate the effective bandwidth in the phase encode direction
BWPerPixPE = BWPerPix / SzOut(2);

% Using the field map, find the new position of each voxel
if (blip_down); Factor = 1; else Factor = -1; end
ynew = y;
xnew = x + Factor * FMap / BWPerPixPE * Sz(2) / SzOut(2);

% Account for aliasing
xnew = mod(xnew,Sz(2));

% Sum/mean in the readout direction (same index for each column)
dx = Sz(2)/SzOut(2);           dy = Sz(1)/SzOut(1);
x_bounds = 0.5:dx:(Sz(2)+0.5); y_bounds = 0.5:dy:(Sz(1)+0.5);
tmp = zeros(SzOut(1),Sz(2)); tmp_x = tmp;
for jj = 1:SzOut(1)
    Idx = ( ynew(:,1) >= y_bounds(jj) ) & ( ynew(:,1) < y_bounds(jj+1) );
    
    if IsSignal % Sum
      tmp(jj,:) = sum(in(Idx,:),1);
    else % Mean
      tmp(jj,:) = mean(in(Idx,:),1);
    end
    
    % Average xnew across this direction also
    tmp_x(jj,:) = mean(xnew(Idx,:),1);
end
in = tmp; xnew = tmp_x;

% Sum/mean in the phase encode direction
for ii = 1:SzOut(2) % Loop over x
  for jj = 1:SzOut(1) % Loop over y
    Idx = ( xnew(jj,:) >= x_bounds(ii) ) & ( xnew(jj,:) < x_bounds(ii+1) );
    
    if IsSignal
      out(jj,ii) = abs(sum(in(jj,Idx),2));
    else
      out(jj,ii) = mean(in(jj,Idx),2);
    end
  end
end

% Unrotate the output
if strcmp(PhaseEncDir,'AP')
   out   = rot90(out,-1);
end

% Taking the mean of no numbers leads to NaN so replace with zeros
out(isnan(out)) = 0;