% This function calculates z-stats, spatial SNR, tSNR, and mean signal
% based on a map of the baseline signal, delta R2*, TE, B0, field map,
% the functional paradigm etc.
%
% Tom Okell, 2011: tokell@fmrib.ox.ac.uk

function [zstat SNR tSNR MeanSig TS Fit] = CalcFMRIStats(opts)

if nargin < 1; opts = []; end

debug = false;

Vars =     {'SigDensity',                                     'dR2star',       'FMap', 'MtxSize', 'SlcThk', 'TE', 'B0', 'Paradigm',    'TR', 'BGMask',                'RegressPhysioNoise', 'HRF',          'HRFt',          'BandWidthPerPix', 'PhaseEncDir', 'blip_down', 'UpSampFactor', 'DropOutMap'  };
Defaults = {CropIm(ones(64,64),[25 25 25 25],true)*10,   ones(64,64)*0.02/30,ones(64,64), 64 ,    5    ,    30 ,   3 ,mod(0:70,20)>9,  3 , [ones(64,1) zeros(64,63)], false ,           load('HRF.txt'), load('HRF_t.txt'),   2000,                'AP',      true       ,  5            , []  };

for ii = 1:length(Vars)
  if isfield(opts, Vars{ii})
     eval([Vars{ii} ' = opts.' Vars{ii} ';'])
  else
     if length(Defaults{ii}) == 1
         str = num2str(Defaults{ii});
     else
         str = 'Too big to display';
     end     
     
     if debug; disp(['Using default value for parameter ' Vars{ii} ': ' str]); end
     
     eval([Vars{ii} ' = Defaults{ii};'])
  end
end

% If no dropout map provided, set to zero
if isempty(DropOutMap)
  DropOutMap = zeros(size(FMap));  
end

% Scale the fieldmap by B0
FMap = FMap * (B0/3);

% Downsample the HRF and normalise
HRFtnew = 0:TR:max(HRFt);
HRF = interp1(HRFt,HRF,HRFtnew);
HRFt = HRFtnew;
HRF = HRF / sum(HRF(:));
if debug; figure; plot(HRFt,HRF); end

% Zero pad so t = 0 is in the centre
% dHRFt = HRFt(2) - HRFt(2);
HRFtnew = -max(HRFt):TR:max(HRFt);
HRF = interp1(HRFt,HRF,HRFtnew);
HRF(isnan(HRF)) = 0; % Replace NaN with 0
HRF = HRF / sum(HRF(:)); % Renormalise
HRFt = HRFtnew;
if debug; hold on; plot(HRFt,HRF,'r'); end

% Convolve the paradigm with the HRF
NVols = length(Paradigm);
t = ((1:NVols) -1)*TR;
EV = conv(double(Paradigm), HRF, 'same');
if debug; figure; plot(t,[Paradigm' EV']); end
%??de-mean here??

% Hardwire FOV and T2* for now
FOV = 220; 
FieldStrengths = [1.5 3 7];
T2stars = [45 30 20];
T2star = interp1(FieldStrengths,T2stars,B0);

% Modify the baseline signal for slice thickness, field
% strength and T2* decay
MeanSig = SigDensity * SlcThk * (B0/3) * exp(-TE/T2star);

% Account for distortion and dropout and resample (also accounts for
% resolution in-plane)
MeanSig = DistortionAndDropout(MeanSig,FMap,UpSampFactor,TE,BandWidthPerPix,MtxSize,SlcThk,PhaseEncDir,blip_down,true ,B0,DropOutMap);
dR2star = DistortionAndDropout(dR2star,FMap,UpSampFactor,TE,BandWidthPerPix,MtxSize,SlcThk,PhaseEncDir,blip_down,false,B0,DropOutMap);

% Calculate the fractional signal change in each voxel, assuming activation
% only exists over 5 mm so increasing beyond this will not improve stats
if SlcThk <= 5;
  Factor = 1;
else
  Factor = 5/SlcThk; 
end
FracSigChange = Factor * TE * dR2star * (B0/3);

% Calculate the pure time series (neglecting the mean signal)
TS = repmat(MeanSig .* FracSigChange,[1 1 NVols]) .* repmat( reshape(EV,1,1,NVols), size(MeanSig));
if debug; dispopt.ref_im = MeanSig; disp_vox(TS,dispopt); end

% Add thermal noise, scaled by bandwidth
ThermalNoiseSD = (MtxSize * BandWidthPerPix) / (64 * 2000) * 10;
TS = TS + randn(size(TS))*ThermalNoiseSD;

% Add physiological noise - cardiac plus random
%?? spatial dependence ??
FieldStrengths = [1.5 3 7];
lambdas = [0.0123 0.0107 0.0086]; % From Triantafyllou, Neuroimage 2005
lambda = interp1(FieldStrengths,lambdas,B0);
PhysiolNoiseSD = lambda * MeanSig;

% Assume 80% of the noise is cardiac related which can be regressed out
CardiacFrac = 0.8;

% Add the "random", non-regressable noise
TS = TS + randn(size(TS)) .* repmat(PhysiolNoiseSD,[1 1 NVols]) * (1-CardiacFrac);

% Add the cardiac noise if not being regressed out
if ~RegressPhysioNoise
  HeartPeriod = 0.95;
  TS = TS + repmat(CardiacFrac*PhysiolNoiseSD,[1 1 NVols]) .* ...
      repmat( reshape(sin(2*pi*t/HeartPeriod),[1 1 NVols] ), size(MeanSig) ) ;
end

if debug; disp_vox(TS,dispopt); end

% ??Smooth data??

% Create mask for running calculations
Mask = MeanSig > 5 * ThermalNoiseSD;

% Convert measured data to matrix (time by voxel)
y = vols2matrix(permute(TS,[1 2 4 3]),Mask)';

% Set up the design matrix, x
x = EV';

% Perform the GLM regression (Modified from code by Mark Jenkinson)
beta = pinv(x)*y;
r = y-x*beta; % Residuals
dof = NVols - rank(x); 
sigma_sq = sum(r.^2,1)/dof; 
var = sigma_sq * inv(x'*x); 
tstat = beta./sqrt(var); 
zstat = t_to_z(tstat,dof);   


% Output zstats and model fit
zstat = matrix2vols(zstat',Mask);
zstat(zstat ==  Inf) =  9; % Can't seem to handle zstats above 9 so fix here
zstat(zstat == -Inf) = -9; % Can't seem to handle zstats above 9 so fix here

Fit = matrix2vols((x*beta)',Mask);
r   = matrix2vols(r',Mask);
% Calculate SNR
% Vol1 = TS(:,:,1);
% BGSD = std(Vol1(Mask));
% SNR = MeanSig / BGSD * 0.65;
SNR = MeanSig / ThermalNoiseSD;

% Calculate tSNR on the residuals
tSNR = MeanSig ./ std(r,[],4);

