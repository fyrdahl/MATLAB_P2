%%
clear
load('smap')

% Pick a set of 16 coils that give a "nice" sensitivity map
coil_set = [2:12 19 23 28 31 32];

smap = smap(:,:,coil_set);

%%
ncoils = size(smap,3);
grappa_factor = 5;

[object] = fast_object;

% Get the object into sensitivity map space
[X,Y] = meshgrid(linspace(0,1,size(object,2)),linspace(0,1,size(object,1)));
[XI,YI] = meshgrid(linspace(0,1,size(smap,2)),linspace(0,1,size(smap,1)));
object = interp2(X,Y,object,XI,YI);

% The data from each coil
image = smap.*repmat(object,[1 1 size(smap,3)]);

% Get the noise correlation matrix and add this to the acquired data
sigma = 20;
rho = 0.2;
Sigma = sigma.*(eye(ncoils)+rho.*(ones(ncoils)-eye(ncoils))); 
Noise = (randn(size(image)) +1i.*randn(size(image)));
[V,D] = eig(Sigma);
W = V*sqrt(D);
for ii=1:size(object,1)
   	for jj=1:size(object,2)
		Noise(ii,jj,:)=W*reshape(Noise(ii,jj,:),[ncoils 1]);
   	end
end

% The acquired data for each coil with noise
image = image+Noise;


% The sampled kspace
full_kspace = ifft2d(image);
undersampled_kspace = zeros(size(full_kspace));
undersampled_kspace(:,1:grappa_factor:end,:) = full_kspace(:,1:grappa_factor:end,:);
sampled_image = fft2d(undersampled_kspace);


% SENSE Reconstruction
full_image = zeros(size(object,1), size(object,2));
for xx = 1:size(sampled_image,1);
     for yy = 1:round(size(sampled_image,2)./grappa_factor);
         
         I = reshape(sampled_image(xx,yy,:),[ncoils 1]);
         PE_points = yy:round(size(object,2)./grappa_factor):size(object,2);
         C = reshape(smap(xx,PE_points,:),[ length(PE_points) ncoils]).';
         full_image(xx,PE_points) = (C'*C)\(C'*I);

     end    
end

imagesc(abs(full_image))
axis image