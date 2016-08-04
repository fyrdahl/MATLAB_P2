close all; clear all;

% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in raw data and display

% get raw data
load brain;
whos

% show raw data
im_data = fft2c(k_data);
whos
show_pair(k_data,im_data);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% effect of distance between k-space lines (delta-k)

% zero every-other line
index = 2:2:256;
k_evenlines = zeros(256,256);
k_evenlines(:,index) = k_data(:,index);
im_evenlines = fft2c(k_evenlines);
show_pair(k_evenlines,im_evenlines);

% keep only even lines
k_evenlinesonly = k_data(:,index);
im_evenlinesonly = fft2c(k_evenlinesonly);
show_pair(k_evenlinesonly,im_evenlinesonly);

% exercise 1: keep only every 3rd line
index = 1:3:256;
k_ex1 = zeros(256,256);
k_ex1(:,index) = k_data(:,index);
im_ex1 = fft2c(k_ex1);
show_pair(k_ex1,im_ex1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% effect of extent of k-space coverage (kmax)

% zero outer 1/2 of k-space
index = 64:192;
k_inner = zeros(256,256);
k_inner(index,index) = k_data(index,index);
im_inner = fft2c(k_inner);
show_pair(k_inner,im_inner);
show_pair(im_data,im_inner);

% keep inner half only
k_inneronly = k_data(index,index);
im_inneronly = fft2c(k_inneronly);
show_pair(k_inneronly,im_inneronly);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partial k-space

% zero upper 5/8 of k-space
index = 1:(256*5/8);
k_partial = zeros(256,256);
k_partial(:,index) = k_data(:,index);
im_partial = fft2c(k_partial);
show_pair(k_partial,im_partial);
show_pair(im_data,im_partial);

% reconstruct with partial k-space filling
index_blank=(256*5/8)+1:256;
index_fill=(256*3/8):-1:1;
k_pk_recon=k_partial;
k_pk_recon(:,index_blank)=real(k_partial(:,index_fill)) ...
                          - i*imag(k_partial(:,index_fill));
im_pk_recon=fft2c(k_pk_recon);
show_pair(im_partial,im_pk_recon);

% reconstruct with homodyne filter
hd_filter=homodyne(256,5/8);
show_pair(k_partial,hd_filter);
k_pk_hd=k_partial.*hd_filter;
im_pk_hd=fft2c(k_pk_hd);
show_pair(im_partial,im_pk_hd);
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPI trajectory

% EPI
load epi_traj;
% L-R, B-T
show_traj(Gx,zeros(size(Gy)));
show_traj(zeros(size(Gx)),Gy);
show_traj(Gx,Gy);
% R-L, B-T
figure; show_traj(-Gx,Gy);
% L-R, T-B
figure; show_traj(Gx,-Gy);
close all;

% spiral
load sp_traj;
show_traj(Gx,Gy);
% flip R/L
figure; show_traj(-Gx,Gy);
close all;

% double-FOV, half-resolution
load epi_traj;
show_traj(Gx,Gy/2);

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')
