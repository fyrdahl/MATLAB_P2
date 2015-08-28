% Fit inversion recovery curve to IR data to determine proton density (PD) and spin-lattice relaxation rate (R1) maps.
%
% FUNCTION [pd r1 eff res] = qmap_t1_fit_ir(data, ti, opts)
%
% Function to fit T1 curve using any IR prepped sequence (ex: IR-EPI, IR-FSE, or IR-SSFSE)
%
% Inputs:
%   data [Np x Nti] - IR-Prepared Magnitude Data (Np: Number of voxels in image,  Nti: Number of Inversion Times)
%   ti   [1  x Nti] - Inversion times (in seconds)
%   opts            - a structure containing optional settings
%       points  - a subset of TI values to use for fitting
%       usemin  - 1 - use all TI values  0 -> Throw out TI value with least signal (null point)
%       fiteff  - 1 - fit additional inversion efficiency 0 - assume 1
%       debug   - debug by plotting data pts and fit curve
%
% Outputs:
%   pd      - proton density factor
%   r1      - Spin-lattice relaxation rate (1/T1)
%   eff     - Inversion pulse efficiency
%   res     - Sum-of-squares residual
%
% Samuel A. Hurley
% University of Wisconsin
% QMAP v2.0 10-Jan-2012
%
% Changelog:
%      v1.0 - initial version based on t2_fit_signa (Mar-2010)
%      v1.1 - fixed initial guess for T1 (1/R1)
%             DO NOT use numerical jacobian for levmar, works faster and
%             better results with analytical version. (Mar-2010)
%      v1.2 - general cleanup so analysis is consistant with t1_fit_ir_vnmr (Mar-2010)
%      v1.3 - added option to NOT fit for inversion efficiency (Apr-2010)
%      v1.4 - Now supports 3-d (mulislice images) due to updated load_dicom
%             command (Jun-2010)
%      v1.5 - check for same R1/R2/TG in all scans (Jun-2010)
%      v2.0 - QMAP Version (Jan-2012)
%__________________________________________________________________________
% Copyright (c) 2005-2012, The Regents of the University of Wisconsin
% All rights reserved. Type 'doc qmap' for full license information.

function [pd r1 eff res] = qmap_t1_fit_ir(data, tivals, opts)

tic;

%% OO. Optomization settings
  % Levenburg-Marquardt /w Numerical Jacobian Matrix, Tolerance is 1e-4
  optim = optimset('Algorithm', 'levenberg-marquardt', 'Jacobian', 'off', 'Tolfun', 1e-4, 'Display', 'off');

%% O. Check input arguments
switch nargin
  case 2
    opts   = struct();    % Use default options
  case 3
    if ~isstruct(opts)
      error('Opts must be a structure.  Type ''help t1_fit_ir_vnmr'' for more information');
    end
  otherwise
    error('You must supply 2 or 3 input arguments.  Type ''help qmap_t1_fit_ir'' for more information');
end

% Grab Options or Set Defaults
if isfield(opts, 'points')
  points = opts.points;
end

if isfield(opts, 'usemin')
  usemin = opts.usemin;
else
  % Default behavior is to use the null point as well
  usemin = 1;
end

if isfield(opts, 'fiteff')
  fiteff = opts.fiteff;
else
  fiteff = 1;
end

if isfield(opts, 'debug')
  debug = opts.debug;
  if debug == 1
    dbgtime = toc;
  end
else
  % Default is to turn off debugging
  debug = 0;
end

% Constants
LOG2 = log(2);

% Banner
disp('==============IR Fit Signa===================')

%% I. Load in image & DICOM header

% Determine which TI points to use
if ~exist('points', 'var')
  % If not specified, use all points
  points = 1:length(tivals);
else
  % Use specified points
  if length(points) > length(tivals)
    error('You specifed more fitting points than avalible TE values');
  end
  
  tivals = tivals(points);
  disp('NOTE: Subset of TI times used for fitting.');
end

disp(['TI Values: ' num2str(tivals,'%01.3f ')]);

% Grab only the TI values specified
data = data(:,points);

% Preallocate outputs
pd  = zeros([size(data,1) 1]);
r1  = zeros([size(data,1) 1]);
eff = zeros([size(data,1) 1]);
res = zeros([size(data,1) 1]);

% Number of non-zero data points
npts = size(find(~(sum(data, 2) == 0)),1);
pt   = 0;
fprintf('T1 IR Fitting:');

%% II. For each point, perform IR fitting
for ii = find(~(sum(data, 2) == 0))';
  % Grab voxel data
  vox_data = data(ii,:);
  
  % Update Progress
  qmap_progressbar(pt/npts);
  pt = pt + 1;
  
  % Initial Guess
  [v idx]   = min(abs(vox_data)); %#ok<ASGLU>
  min_ti    = tivals(idx);
  x         = [max(vox_data(:)) LOG2/min_ti 1];  % TInull = T1/ln(2)
  
  % Check if we should use the minimum-singal TI point
  if usemin == 0
    % Remove the TI value and corresponding data point
    ti = tivals;
    ti(idx) = [];
    vox_data(idx) = [];
  else
    % Use all TI values
    ti = tivals;
  end
  
  % Minimize /w Levenburg-Marquardt
  [x residual] = lsqnonlin(@ir_model, x, [], [], optim);
  
  % Save results
  pd(ii)       = x(1);
  r1(ii)       = x(2);
  if fiteff == 1
    eff(ii)    = x(3);
  else
    eff(ii)    = 1;
  end
  res(ii)      = residual;
  
  % Optinal debug plotting
  if debug == 1 && (toc - dbgtime) > .250
    % Only re-draw every 1 second
    
    dbgtime = toc;

    % Setup Figure
    if ~exist('dbfig', 'var')
      dbfig = figure;
    end
 
    % Setup vector
    ti = min(tivals):.01:max(tivals);

    % Plot actual vs fitted data
    figure(dbfig);
    plot(tivals, data(ii,:), 'o', ti, ir_model(x,1), '-');

    % Figure captions & Legend
    title(['T1: ' num2str(1/x(2), '%0.2f') ' s']);
    legend('MRI Data', 'Fit Curve');
    xlabel('Inversion Time [s]');
    ylabel('Signal [a.u.]');
    drawnow;
  end
  
end

% Done.
qmap_progressbar(1);
toc;

%% V. IR Model Function
  function [F J] = ir_model(x, dbgmode) %#ok<INUSD>
    % Pull out parameters
    pd_mod    = x(1);
    r1_mod    = x(2);
    
    if fiteff == 1
      eff_mod = x(3);
    else
      % Don't fit inversion efficiency
      eff_mod = 1;
    end
    
    % MR Signal Mz (z-magnetization)
    mz = pd_mod.*(1-2.*eff_mod.*exp(-ti.*r1_mod));
    
    % Output vector
    F = abs(mz);      % Use sign-sensitive signal
    
    if ~exist('dbgmode', 'var')
      F = F - vox_data; % Difference vector
    end
    
    % Jacobian / First Derivatives of Mz
    if nargout == 2
      dmz_dm0 = 1 - 2.*eff_mod.*exp(-ti.*r1_mod);
      dmz_dr1 = 2.*pd_mod.*eff_mod.*ti.*exp(-ti.*r1_mod);
      dmz_def = -2.*pd_mod.*exp(ti.*r1_mod);
      
      % Output jacobian
      J = [(dmz_dm0)' (dmz_dr1)' (dmz_def)'];
    end
  end

end