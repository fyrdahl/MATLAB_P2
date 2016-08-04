
ig = image_geom('nx', 32, 'ny', 28, 'dx', 1, 'offsets', 'dsp');
ig.mask = ig.circ(1+ig.nx/2, 1+ig.ny/2) > 0;
N = ig.dim;
J = [6 6];
nufft_args = {N, J, 2*N, N/2, 'table', 2^10, 'minmax:kb'};

Gn = Gnufft(ig.mask, {omega, nufft_args{:}});
Gm = Gmri(kspace, ig.mask, ...
    'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);

Tm = build_gram(Gm, 1);
x1 = ig.embed(Gm' * (Gm * x0(ig.mask)));
x2 = ig.embed(Tm * x0(ig.mask));
im clf, im(stackup(x1,x2))


[oo1 oo2] = ndgrid(	2*pi*([0:N(1)-1]/N(1) - 0.5), ...
2*pi*([0:N(2)-1]/N(2) - 0.5));
yd_g = griddata(omega(:,1), omega(:,2), yi, oo1, oo2, 'cubic');
yd_g(isnan(yd_g)) = 0;

xg = ifft2(fftshift(yd_g));
im(4, abs(xg), '|x| "gridding"', clim), cbar








