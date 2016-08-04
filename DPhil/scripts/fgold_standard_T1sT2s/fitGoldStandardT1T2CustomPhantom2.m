%% Gold Standard T1 and T2 maps of SphereD170

%% 1. Set up paths
%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/DPhil';
workingdir = '/Users/jallen/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze
%
savingdir = '/Users/jallen/Documents/DPhil';
addpath(genpath(savingdir));

% If working on jalapeno00, uncomment the following lines:
% addpath(genpath('/Applications/fsl/'))
% addpath(genpath('/usr/local/fsl/bin'))
% addpath(genpath('/opt/fmrib/fsl/etc/matlab'))

%% 2. Provide Scan Details and Read Data
ID = 'MR_PHYSICS_439';
date = '20160311';
TIprotocol = 'ep2d_t1_TI';
TEprotocol = 'ep2d_t2_TE';
TI =  [100:50:600,700:100:1100];
TE = [25:10:175];
TIseriesStart = 2;
TEseriesStart = 18;
%
for n = 1:6
TIdata(:,:,n) = nifti2mat(['/Users/jallen/Documents/DPhil/data/NIfTI/20160311_MR_PHYSICS_439/TIgroup/',date,'_00',num2str(TIseriesStart+(n-1)),'_',TIprotocol,num2str(TI(n)),'.nii']);
end
for n = 9:numel(TI)
TIdata(:,:,n) = nifti2mat(['/Users/jallen/Documents/DPhil/data/NIfTI/20160311_MR_PHYSICS_439/TIgroup/',date,'_0',num2str(TIseriesStart+(n-1)),'_',TIprotocol,num2str(TI(n)),'.nii']);
end
for n = 1:numel(TE)
TEdata(:,:,n) = nifti2mat(['/Users/jallen/Documents/DPhil/data/NIfTI/20160311_MR_PHYSICS_439/TEgroup/',date,'_0',num2str(TEseriesStart+(n-1)),'_',TEprotocol,num2str(TE(n)),'.nii']);
end

%% 3. Fit T1 curves to Inversion recovery data
TIdata = reshape(TIdata,[size(TIdata,1)*size(TIdata,2),size(TIdata,3)]);
opts = struct('debug',0,'fiteff',1);
disp('Fitting T1...')
for n = 1:size(TIdata,1)
[T1fitpd(n) r1(n) T1fitBeta(n)] = qmap_t1_fit_ir(TIdata(n,:), 0.001*TI,opts); % Sam Hurley's function
    T1fitmap(n,1) = (1/r1(n))*1000;
       if mod(100*n/size(TIdata,1),1) == 0;
    disp(['Fitting T1...progress:',num2str(100*n/size(TIdata,1)),'%'])
       end
end
T1fitmap = reshape(T1fitmap,64,64);
T1fitBeta = reshape(T1fitBeta,64,64);
T1fitpd = reshape(T1fitpd,64,64);
disp('Fitting T1...Finished.')

%% 4. Fit T2 curves to spin echo images
disp('Fitting T2...')
TEdata = reshape(TEdata,[size(TEdata,1)*size(TEdata,2), size(TEdata,3)]);
T2fitmap = zeros(1,size(TEdata,1));
T2model = @(a,T2,c, x) a*exp(-x*(1/T2)) + c;
for n = 1:size(TEdata,1)
    [T2Curve] = fit(TE',TEdata(n,:)',T2model,'Upper',[5000 5000 2500],'Lower',[0 0 0],'StartPoint',[1000 100 100]);
    T2fitmap(n,1) = T2Curve.T2;
    T2fita(n) = T2Curve.a;
    T2fitc(n) = T2Curve.c;
    if mod(100*n/size(TEdata,1),1) == 0;
    disp(['Fitting T2...progress:',num2str(100*n/size(TEdata,1)),'%'])
    end
end
T2fitmap = reshape(T2fitmap,64,64);
T2fita = reshape(T2fita,64,64);
T2fitc = reshape(T2fitc,64,64);
disp('Fitting T2...Finished.')

%% 5. Plot Results
%
figure, 
subplot 231
imagesc(T1fitmap)
c = colorbar;
caxis([0 2000])
%
subplot 232
imagesc(T1fitBeta)
%
subplot 233
imagesc(T1fitpd)
%
subplot 234
imagesc(T2fitmap)
c = colorbar;
caxis([0 1500])
%
subplot 235
imagesc(T2fita)
%
subplot 236
imagesc(T2fitc)


%% 6. Save Results
dt = datetime('now');
save([savingdir,'/MAT-files/T1fitmap',datestr(dt),'.mat'],'T1fitmap')
save([savingdir,'/MAT-files/T1fitBeta',datestr(dt),'.mat'],'T1fitBeta')
save([savingdir,'/MAT-files/T1fitpd',datestr(dt),'.mat'],'T1fitpd')
save([savingdir,'/MAT-files/T2fitmap',datestr(dt),'.mat'],'T2fitmap')
save([savingdir,'/MAT-files/T2fita',datestr(dt),'.mat'],'T2fita')
save([savingdir,'/MAT-files/T2fitc',datestr(dt),'.mat'],'T2fitc')
