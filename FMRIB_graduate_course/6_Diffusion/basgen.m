function s=basgen(b,v,s0,d,f,th,ph,dstd)
% s=pvm_forwardmodel(bvals,bvecs,s0,d,f,th,ph,dstd)
%
% bvals is 1xN
% bvecs is 3xN
% s0, d are Px1  (P=number of voxels)
% f,th,ph are PxF (F=number of fibres)
%
% Output is PxN
if(nargin<8)
    model=1;
else
    model=2;
    dstd=repmat(dstd,1,size(f,2));
    valid=dstd(:,1)<1e-5;

end
P=size(s0,1);
N=size(b,2);
s=zeros(P,N);
[fx,fy,fz]=sph2cart(ph,pi/2-th,1);
d=repmat(d,1,size(f,2));
for i=1:N        
    bg=b(i)*(fx*v(1,i)+fy*v(2,i)+fz*v(3,i)).^2;    
    if(model==1)
        s(:,i)=s0.*(sum( f.*exp(-bg.*d) ,2)+(1-sum(f,2)).*exp(-b(i)*d(:,1)));
    else        
        alph=(d.^2)./(dstd.^2);
        beta=d./(dstd.^2);
        s(:,i)=s0.*(sum(f.*exp(alph.*log(beta./(beta+bg))),2) +(1-sum(f,2)).*exp(alph(:,1).*log(beta(:,1)./(beta(:,1)+b(i)))));
        
        s(~valid,i)=s0(~valid).*(sum( f(~valid,:).*exp(-bg(~valid,:).*d(~valid,:)) ,2)+(1-sum(f(~valid,:),2)).*exp(-b(i)*d(~valid,1)));
    end
end
