function res = load_bedpostx_results(bpxdir)
% res = load_bedpostx_results(bpxdir)
% 
% addpath([getenv('FSLDIR') '/etc/matlab']);
nfibres=1;
while(exist([bpxdir '/dyads' num2str(nfibres) '.nii.gz'],'file'))
    nfibres=nfibres+1;
end
nfibres=nfibres-1;

% res.mask=read_avw([bpxdir '/mean_f1samples'])>0;
mm = load_nii([bpxdir '/mean_f1samples.nii.gz']);
res.mask = mm.img > 0;



N=numel(res.mask);

% res.s0=reshape(read_avw([bpxdir '/mean_S0samples']),N,1);
mm = load_nii([bpxdir '/mean_S0samples.nii.gz']);
mm = mm.img;
res.s0=reshape(mm,N,1);

% res.d=reshape(read_avw([bpxdir '/mean_dsamples']),N,1);
mm = load_nii([bpxdir '/mean_dsamples.nii.gz']);
mm = mm.img;
res.d=reshape(mm,N,1);



if(exist([bpxdir '/mean_d_stdsamples.nii.gz'],'file'))
%     res.dstd=reshape(read_avw([bpxdir '/mean_d_stdsamples']),N,1);
    mm = load_nii([bpxdir '/mean_d_stdsamples.nii.gz']);
    mm = mm.img;
    res.dstd = reshape(mm,N,1);
    
end
for i=1:nfibres
    
%     res.f(:,i)=reshape(read_avw([bpxdir '/mean_f' num2str(i) 'samples']),N,1);
    mm = load_nii([bpxdir '/mean_f' num2str(i) 'samples.nii.gz']);
    mm = mm.img;
    res.f(:,i)=reshape(mm,N,1);
    
%     dyads=reshape(read_avw([bpxdir '/dyads' num2str(i)]),N,3);
    mm = load_nii([bpxdir '/dyads' num2str(i) '.nii.gz']);
    mm = mm.img;
    dyads=reshape(mm,N,3);
    
    
    x=dyads(:,1);y=dyads(:,2);z=dyads(:,3);
    [ph,th]=cart2sph(x,y,z);th=pi/2-th;
    res.th(:,i)=th;
    res.ph(:,i)=ph;    
end
