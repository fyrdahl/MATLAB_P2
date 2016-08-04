function res = fast_dtifit_eig(y,b,v,y0)
% res = fast_dtifit(data,bvals,bvecs)
% data is PxN
% bvals is 1xN
% bvecs is 3xN
%
% S.Jbabdi 09/14
% mod by wenchuan

[nx,ny,N] = size(y);
P = nx * ny;
% 
% [P,N]=size(y);

y = reshape(y,[P,N]);

y0 = reshape(y0,[P,N]); % use good img for mask
mask=sum(abs(y0(:,1)),2)>1e-3;
y=y(mask,:);

res.pred = zeros(P,N);
res.S0   = zeros(P,1);
res.L1   = zeros(P,1);
res.L2   = zeros(P,1);
res.L3   = zeros(P,1);
res.MD   = zeros(P,1);
res.FA   = zeros(P,1);
res.V1   = zeros(P,3);

v=v';b=b';
A = [v(:,1).*v(:,1).*b ...
    2*v(:,1).*v(:,2).*b ...
    2*v(:,1).*v(:,3).*b ...
    v(:,2).*v(:,2).*b ...
    2*v(:,2).*v(:,3).*b ...
    v(:,3).*v(:,3).*b ...
    ones(length(b),1)];
r = -pinv(A)*log(y'); % r is 7xN

res.pred(mask,:) = exp(-A*r)';
res.S0(mask)=exp(-r(7,:))';


% explicit roots
a00 = r(1,:)';
a01 = r(2,:)';
a02 = r(3,:)';
a11 = r(4,:)';
a12 = r(5,:)';
a22 = r(6,:)';

% eig is very fast

tensor = zeros(size(y,1),3,3);
eigvals = zeros(size(tensor));
eigvecs = zeros(size(tensor));

tensor(:,1,1) = a00;
tensor(:,1,2) = a01;
tensor(:,1,3) = a02;
tensor(:,2,1) = a01;
tensor(:,2,2) = a11;
tensor(:,2,3) = a12;
tensor(:,3,1) = a02;
tensor(:,3,2) = a12;
tensor(:,3,3) = a22;

% wenchuan add, no idea about the consequence
tensor(find(isnan(tensor))) = 0;
tensor(find(isinf(tensor))) = 0;

for ii = 1 : size(y,1)
   
    [eigvecs(ii,:,:), eigvals(ii,:,:)] = eig(squeeze(tensor(ii,:,:)));    
    
end

root0 = squeeze(eigvals(:,1,1));
root1 = squeeze(eigvals(:,2,2));
root2 = squeeze(eigvals(:,3,3));

roots=sort([root0 root1 root2],2,'descend');
res.L1(mask)=roots(:,1);
res.L2(mask)=roots(:,2);
res.L3(mask)=roots(:,3);
res.MD(mask)=mean(roots,2);
res.FA(mask)=sqrt(3/2).*sqrt(sum((roots-repmat(res.MD(mask),1,3)).^2,2))./sqrt(sum(roots.^2,2));

res.V1(mask,:) = squeeze(eigvecs(:,:,3));

% reshape
res.MD = reshape(res.MD,[nx,ny]);

% found that the mask based on dwi data is not enough to select background,
% so use MD
msk = zeros(size(res.MD));
msk(find(res.MD > 1e-4)) = 1;
res.FA = reshape(res.FA,[nx,ny]);

res.FA = res.FA .* msk;

res.V1 = reshape(res.V1,[nx,ny,3]);
res.V1 = res.V1 .* repmat(msk,[1,1,3]);









% 
% V = reshape(res.V1,[145,174,3]);
% 
% va = abs(V);
% % image(V/max(V(:)))
% image(va/max(va(:)))
% FA = reshape(res.FA,[145,174]);
% va = abs(V.*repmat(FA,[1,1,3]));
% image(va/max(va(:)))

 
