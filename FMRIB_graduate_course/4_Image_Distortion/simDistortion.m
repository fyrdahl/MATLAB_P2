% Distortion and chemical shift simulator for the FMRIB grad course
%
%   jmeakin@fmrib.ox.ac.uk - November 2011
%
% Options, you can set as many of these as you like!
% Arguments that take a number:
%     GRAPPA - set the GRAPPA factor (integer, ie. 2 or 3), default 1
%     Bandwidth - set the bandwidth in Hz/pixel, default 2000
%     FieldStrength -  set the field sterngth in T, default 3
%     Matrix - set the matrix size, default 64
%     SliceThickness - set the slice thickness in mm, default 3
%     TE - set the TE in s, detault 30
%
% Arguments that take a bool:
%     AddNoise - add in some rician noise, default true
%     FatSat - perform fat sat, default false
%     BlipUp - PE blips up, default false
%     Dropout - simulate dropout, default true
%     AddSinus - add in field offsets due to sinuses, default false
%     PhaseEncodeLR - phase encode left right, default false
%
% Useage:
%
%     simDistortion();
%         Runs the simulator with the default parameters
%
%     options.Bandwidth = 1000;
%     simDistortion(options);
%         Runs the simulation with a bandwidth of 1000 Hz/px
%
%     options.Bandwidth = 1000;
%     options.FatSat = true;
%     options.PhaseEncodeLR = true;
%     simDistortion(options);
%         Runs the simulation with a bandwidth of 1000 Hz/px, fat sat on
%         with phase encoding left-right
%



function [out] = simDistortion(options)

% if (nargin>0)
%   S = fieldnames(orderfields(options));
%   s = lower(S);
%  [A,idx1] = unique(s,'first');
%  [~,idx2] = unique(s,'last');
%  if any(idx1~=idx2)
%     opt1=S(idx1(find(idx2~=idx1)));
%     opt2=S(idx2(find(idx1~=idx2)));
%     error(['Multiple option definitions ' cell2mat(opt1) ' and ' cell2mat(opt2) ...
%            '. Please remove one definition by using:  ' ...
%            ' options=rmfield(options,''' cell2mat(opt1) ''')'])
%  end
% end

if (nargin<1) options=default_options();
else          options=default_options(options);
end;

if (~isnumeric(options.Bandwidth)); error('The value of options.Bandwidth must be a number.'); end;
LR_BW_px = options.Bandwidth;
if (~isnumeric(options.GRAPPA)); error('The value of options.GRAPPA must be a number.'); end;
% KM - faking GRAPPA for now, just reducing bandwidth
grappa_fact = options.GRAPPA;
%LR_BW_px = LR_BW_px*options.GRAPPA;
%grappa_fact = 1;
if (~isnumeric(options.FieldStrength)); error('The value of options.FieldStrength must be a number.'); end;
field = options.FieldStrength;
if (~isnumeric(options.Matrix)); error('The value of options.Matrix must be a number.'); end;
matrix = options.Matrix;
if (~isnumeric(options.SliceThickness)); error('The value of options.SliceThickness must be a number.'); end;
thk = options.SliceThickness * 1e-3;
if (~isnumeric(options.TE)); error('The value of options.TE must be a number.'); end;
TE = options.TE * 1e-3;

if (~length(find([0 1])==options.AddNoise)); error('The value of options.AddNoise must be logical, use true or false.'); end;
add_noise = options.AddNoise;
if (~length(find([0 1])==options.FatSat)); error('The value of options.FatSat must be logical, use true or false.'); end;
fat_sat = options.FatSat;
if (~length(find([0 1])==options.BlipUp)); error('The value of options.BlipUp must be logical, use true or false.'); end;
blip_down = ~options.BlipUp;
if (~length(find([0 1])==options.Dropout)); error('The value of options.Dropout must be logical, use true or false.'); end;
dropout = options.Dropout;
if (~length(find([0 1])==options.AddSinus)); error('The value of options.AddSinus must be logical, use true or false.'); end;
sim_sinus = options.AddSinus;
if (~length(find([0 1])==options.PhaseEncodeLR)); error('The value of options.PhaseEncodeLR must be logical, use true or false.'); end;
rotate_PE = options.PhaseEncodeLR;


%% ~~~~~~~~~~~~ Defaults
FOV = 0.3; % Phantom size, m
RO_FOV = FOV; % Readout FOV, we'll never undersample this
gamma = 2.6752220 * 10 ^ 8; % gyromagnetic ratio
max_field_map_shift = 100; % Hz/T

% Change the size of the phantom so we can perform grappa
%if ~(mod(matrix,grappa_fact) == 0)
%    matrix = ceil(matrix/grappa_fact).*grappa_fact;
%end

T_lin = 1./(LR_BW_px); % time to acquire one line, s

G = 2.*pi.*LR_BW_px.*matrix./(gamma.*RO_FOV); % Check RO gradient strength, T/m
%if(G > 40e-3) % System gradient maximum, T/m
if(G/options.GRAPPA > 40e-3) % System gradient maximum, T/m
    error('Our system isn''capable of producing a Gradient > 40 mT/m, try reducing the Bandwidth or reducing Resolution (FOV/Matrix)')
end

Slew_rate_max = 170; % T/m/s
rise_time = G./Slew_rate_max; % s

T_esp = (T_lin + 2.*rise_time)./grappa_fact; % Echo spacing, s

n_PE_lines = matrix; % Number of PE lines
T_etl = n_PE_lines.*T_esp; % Total echo train time, s

% NOTE: This is not the case for partial fourier...
AP_BW_px = 1./T_etl; % Bandwidth in PE direction

% Set the LR and AP Bandwidths (Hz)
LR_BW = LR_BW_px .* matrix;
AP_BW = AP_BW_px .* n_PE_lines;

% We'll bin the object into Subsample bins in the PE direction
Subsample = 10;
%object = fast_phantom(Subsample .* matrix);
%fat_mask = object >= 0.5;
[object, fat_mask] = fast_object;

% Make the object square
padding_size = round(abs(size(object) - max(size(object)))/2);
object = padarray(object,padding_size);
fat_mask = padarray(fat_mask,padding_size);

% Upscale the object
object = abs(fft2d(padarray(ifft2d(object),((Subsample*matrix)-size(object))./2.0)));
object = object./max(object(:));
fat_mask = abs(fft2d(padarray(ifft2d(fat_mask),((Subsample*matrix)-size(fat_mask))./2.0)));
fat_mask = fat_mask > 0.1;

if (rotate_PE)
    object = rot90(object);
    fat_mask = rot90(fat_mask);
end

% Set the PE Field of view, n_PE_lines
PE_FOV = RO_FOV ;% ./ grappa_fact ;

% Remove the fat
if (fat_sat)
    object(fat_mask) = 0.1*object(fat_mask);
end
%%

% Make a matrix of absolute positions (m)
abs_pos = zeros(Subsample .* matrix, Subsample .* matrix,2);
abs_pos(:,:,1) = repmat(linspace(0,FOV, Subsample .* matrix) , Subsample .* matrix,1);
abs_pos(:,:,2) = repmat(linspace(0,FOV, Subsample .* matrix)',1, Subsample .* matrix);

% Create a field map (Hz)
field_map = zeros(Subsample .* matrix);
if (sim_sinus)
    % Add in three gaussian functions for the sinuses
    [X,Y] = meshgrid(1:size(field_map,1),1:size(field_map,2));
    field_map = field_map + exp(-((X-0.5.*size(field_map,1)).^2 + (Y-0.15.*size(field_map,2)).^2)./(15*size(field_map,1)));
    field_map = field_map + exp(-((X-0.25.*size(field_map,1)).^2 + (Y-0.5.*size(field_map,2)).^2)./(5*size(field_map,1)));
    field_map = field_map + exp(-((X-0.75.*size(field_map,1)).^2 + (Y-0.5.*size(field_map,2)).^2)./(5*size(field_map,1)));
end

% Spatial distortion? Set to zero for now
field_map =  field_map +...
    0.*(abs_pos(:,:,1)- FOV/2) ...
    + 0.*(abs_pos(:,:,2)- FOV/2) ...
    + 0.*(abs_pos(:,:,1)- FOV/2).*(abs_pos(:,:,2)- FOV/2) ...
    + 0.*((abs_pos(:,:,1) - FOV/2).^2 - (abs_pos(:,:,2)- FOV/2).^2);

field_map = field_map - mean(field_map(:)); % De-mean (shim)
field_map = field .* max_field_map_shift .* field_map;

% Rotation for LR PE encoding
if (rotate_PE)
    field_map = rot90(field_map);
end

% Dropout
if (dropout)
    object = object.*exp(1i.*field_map.*TE);
    % Through plane dephasing
    object = object.*exp(1i.*field_map.*TE.*thk.*5);
end

% Bin in RO
object_new = zeros(Subsample .* matrix,matrix);
field_map_new = zeros(Subsample .* matrix,matrix);
fat_mask_new = false(Subsample .* matrix,matrix);
for ii = 1:matrix
    object_new(:,ii) = sum(object(:,(((ii-1)*Subsample)+1):(ii*Subsample) ),2);
    field_map_new(:,ii) = mean(field_map(:,(((ii-1)*Subsample)+1):(ii*Subsample) ),2);
    fat_mask_new(:,ii) = any(fat_mask(:,(((ii-1)*Subsample)+1):(ii*Subsample) )');
end

% Create the MPRAGE image
MPRAGE = object_new;

% Add the fat to the field map
field_map = field_map_new;
if (~fat_sat)
    fat_offres = field .* 3.5e-6 .* gamma ./ (2*pi);
    field_map(fat_mask_new) = field_map(fat_mask_new) + fat_offres;
end


% Get the positions of the voxel
abs_pos = zeros(Subsample .* matrix,matrix,2);
abs_pos(:,:,1) = repmat(linspace(0,FOV,  matrix) , Subsample .*matrix,1);
abs_pos(:,:,2) = repmat(linspace(0,FOV, Subsample .* matrix)',1,  matrix);

% Using the field map, find the new position of the voxel
new_pos = zeros(Subsample .* matrix,matrix,2);
new_pos(:,:,1) = abs_pos(:,:,1) + field_map .* FOV ./ LR_BW;
if (blip_down)
    new_pos(:,:,2) = abs_pos(:,:,2) + field_map .* PE_FOV ./ AP_BW;
else
    new_pos(:,:,2) = abs_pos(:,:,2) - field_map .* PE_FOV ./ AP_BW;
end

% If we're doing GRAPPA, center the image (looks nicer...)
%if grappa_fact > 1
%    new_pos(:,:,2) = new_pos(:,:,2) - PE_FOV ./ grappa_fact;
%end

% Sort the new positions of the voxel
intensities = zeros(Subsample .* matrix,matrix);
B = zeros(Subsample .* matrix,matrix);
for ii = 1:size(new_pos,2)
    [B(:,ii),IX] = sort(mod(new_pos(:,ii,2),PE_FOV));
    tmp = object_new(:,ii);
    intensities(:,ii) = tmp(IX);
end

% Bin these into our sampling points
sampling = linspace(0,PE_FOV,matrix + 1);
new_sample = zeros(matrix,matrix);
for ii = 1:(length(sampling) - 1)
    for jj = 1:size(new_pos,2);
        idx = (B(:,jj) <= sampling(ii+1)) & (B(:,jj) > sampling(ii));
        new_sample(ii,jj) = abs(sum(intensities(idx,jj)));
    end
end

% Bin the GRE
simple_GRE = zeros(matrix, matrix);
for ii = 1:matrix
    simple_GRE(ii,:) = abs(sum(MPRAGE((((ii-1)*Subsample)+1):(ii*Subsample),: ),1));
end


%%

% Get the sensitivity map
smap = get_smap;
smap = reshape(smap,[64 64 16]);
ncoils = size(smap,3);
new_sample = smap.*repmat(new_sample,[1 1 size(smap,3)]);

full_kspace = ifft2d(new_sample);

if(add_noise)
    vox_vol = (FOV./matrix).^2.*thk;
    % Get the noise correlation matrix and add this to the acquired data
    sigma = 1.0E-11*sqrt(LR_BW_px)./(field.*vox_vol);
    rho = 0.2;
    Sigma = sigma.*(eye(ncoils)+rho.*(ones(ncoils)-eye(ncoils)));
    Noise = (randn(size(full_kspace)) +1i.*randn(size(full_kspace)));
    [V,D] = eig(Sigma);
    W = V*sqrt(D);
    for ii=1:size(full_kspace,1)
        for jj=1:size(full_kspace,2)
            Noise(ii,jj,:)=W*reshape(Noise(ii,jj,:),[ncoils 1]);
        end
    end
    
    % The acquired data for each coil with noise
    full_kspace = full_kspace+Noise;
end

% The sampled kspace
undersampled_kspace = zeros(size(full_kspace));
undersampled_kspace(:,1:grappa_fact:end,:) = full_kspace(:,1:grappa_fact:end,:);
sampled_image = fft2d(undersampled_kspace);

% SENSE Reconstruction
new_sample = zeros(size(new_sample,1), size(new_sample,2));
for xx = 1:size(sampled_image,1);
    for yy = 1:round(size(sampled_image,2)./grappa_fact);
        
        I = reshape(sampled_image(xx,yy,:),[ncoils 1]);
        PE_points = yy:round(size(new_sample,2)./grappa_fact):size(new_sample,2);
        C = reshape(smap(xx,PE_points,:),[ length(PE_points) ncoils]).';
        new_sample(xx,PE_points) = (C'*C)\(C'*I);
        
    end
end
new_sample = abs(new_sample);


% Rotation for LR PE encoding
if (rotate_PE)
    field_map = rot90(field_map,3);
    new_sample = rot90(new_sample,3);
    simple_GRE = rot90(simple_GRE,3);
end


% Plotting

warning('off','MATLAB:HandleGraphics:noJVM')

figure()
%set(gcf,'WindowStyle','docked')
clf
subplot(2,2,1)
imagesc(fast_object)
axis image
axis off
title('Actual Brain')

subplot(2,2,2)
imagesc(field_map)
axis square
axis off
title('Field Map')

colormap gray

%figure(2)
%clf
subplot(2,2,3)
imagesc(simple_GRE,[0 0.9.*Subsample.^2])
axis image
axis off
title('Simple GRE')


subplot(2,2,4)
imagesc(new_sample)
axis image
axis off
title('EPI')

% RT = round(rise_time*1e6);
% TL = round(T_lin*1e6);
% subplot(2,2,4)
% plot([0, RT, RT+TL, 2*RT+TL, 3*RT+TL, 3*RT+2*TL, 4*RT+2*TL],[0 G G 0 -G -G 0].*1e3)
% axis([0 inf -40 40])
% colormap gray

%out.Min_TE = 1e3.*T_etl./2;
out.PE_Spacing = 1e6.*T_esp;
out.Rise_Time = rise_time*1e6;
out.RO_Time = T_lin*1e6;
out.SliceTR = ((T_etl./(2*grappa_fact) + TE) + 10E-3);
%out.EPI = new_sample;



end
