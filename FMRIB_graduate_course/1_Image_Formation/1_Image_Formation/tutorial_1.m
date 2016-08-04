load brain
load epi_traj

%	Appendix A1
%	==============================================================

im	=	ifft2c(k_data);

show_pair(k_data, im)

A	=	zeros(64, 32);

show_traj(Gx, Gy)

%	Part 1 - K-space Sampling
%	==============================================================

%	Brain
show_pair(k_data, ifft2c(k_data));

%	Noise
n 	=	randn(128);
show_pair(fft2c(n), n);

%	Rectangle
show_pair(fft2c(rect), rect);

%	Circle
show_pair(fft2c(circ), circ);

%	Gaussians
show_pair(fft2c(gauss_fat), gauss_fat);
show_pair(fft2c(gauss_thin), gauss_thin);


%	Conjugate Symmetry
disp(k_data(127:131, 127:131));

%	Shift vs. Magnitude
show_pair(fft2c(im), im);
show_pair(fft2c(circshift(im, 128)), circshift(im, 128));

%	Zero-mean image
show_pair(fft2c(im-mean(im(:))), im-mean(im(:)));


%	Part 1.2 - ∆k and FOV
%	==============================================================

%	Sub-sampling ky by 2x
k_sub	=	k_data(:, 2:2:end);
x_sub	=	ifft2c(k_sub);
show_pair(k_sub, x_sub);

%	Sub-sampling ky by 3x
k_sub	=	k_data(:, 3:3:end);
x_sub	=	ifft2c(k_sub);
show_pair(k_sub, x_sub);

%	Part 1.3 - k_max and resolution
%	==============================================================

%	Restrict kx
k_subx	=	k_data(64+(1:128), :);
show_pair(k_subx, ifft2c(k_subx), 1);

%	Restrict ky
k_suby 	= 	k_data(:, 64+(1:128));
show_pair(k_suby, ifft2c(k_suby), 1);

%	Restrict kx & ky
k_subxy = 	k_data(64+(1:128), 64+(1:128));
show_pair(k_subxy, ifft2c(k_subxy), 1);

%	Zero-pad Images
k_xy_zeropad = zeros(256);
k_xy_zeropad(65:192, 65:192)	= 	k_subxy;
show_pair(k_xy_zeropad, ifft2c(k_xy_zeropad));

%	Part 2 - K-space Trajectories
%	==============================================================

%	EPI
show_traj(Gx, Gy);

%	Spiral
Gx_sp 	= 	(1:1000).*cos((pi/100)*(1:1000));
Gy_sp 	= 	(1:1000).*sin((pi/100)*(1:1000));
show_traj(Gx_sp, Gy_sp);

%	Stochastic
Gx_rand	=	randn(1, 1000);
Gy_rand = 	randn(1, 1000);
show_traj(Gx_rand, Gy_rand);

%	Hybrid Spiral-Stochastic
show_traj(Gx_sp/1000, Gy_rand);

%	Rosette
a 	=	17;
b 	=	23;
t 	=	0:0.005:4;
G_rose = (a+b)*exp(1j*(a+b)*t)+(a-b)*exp(-1j*(a-b)*t);
show_traj(real(G_rose), imag(G_rose));