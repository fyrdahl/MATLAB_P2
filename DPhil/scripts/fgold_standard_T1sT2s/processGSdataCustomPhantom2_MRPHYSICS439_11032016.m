%% Gold Standard T1 and T2 maps of SphereD170

%% 1. Set up paths
%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/DPhil';
workingdir = '/Users/jallen/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze


savingdir = '/Users/jallen/Documents/DPhil';
addpath(genpath(savingdir));

% If working on jalapeno00, uncomment the following lines:
% addpath(genpath('/Applications/fsl/'))
% addpath(genpath('/usr/local/fsl/bin'))
% addpath(genpath('/opt/fmrib/fsl/etc/matlab'))

%% 2. Load images

% Scan Details
ID = 'MR_PHYSICS_439';
date = '20160311';
TIprotocol = 'ep2d_t1_TI';
TEprotocol = 'ep2d_t2_TE';
TI =  [100:50:600,700:100:1100];
TE = [25:10:175];
TIseriesStart = 2;
TEseriesStart = 18;

for n = 1:6
    TIdata(:,:,n) = read_avw(['/Users/jallen/Documents/DPhil/data/NIfTI/20160311_MR_PHYSICS_439/TIgroup/',date,'_00',num2str(TIseriesStart+(n-1)),'_',TIprotocol,num2str(TI(n)),'.nii']);
end
for n = 9:numel(TI)
    TIdata(:,:,n) = read_avw(['/Users/jallen/Documents/DPhil/data/NIfTI/20160311_MR_PHYSICS_439/TIgroup/',date,'_0',num2str(TIseriesStart+(n-1)),'_',TIprotocol,num2str(TI(n)),'.nii']);
end

for n = 1:numel(TE)
    TEdata(:,:,n) = read_avw(['/Users/jallen/Documents/DPhil/data/NIfTI/20160311_MR_PHYSICS_439/TEgroup/',date,'_0',num2str(TEseriesStart+(n-1)),'_',TEprotocol,num2str(TE(n)),'.nii']);
end

%% 3. Fit T1 and T2 curves to Gold Standard images
TIdata = reshape(TIdata,[size(TIdata,1)*size(TIdata,2),size(TIdata,3)]);
opts = struct('debug',0,'fiteff',1);

disp('Fitting T1...')
for n = 1:size(TIdata,1)
    [pd(n) r1(n) eff(n)] = qmap_t1_fit_ir(TIdata(n,:), 0.001*TI,opts);
    BetaMap(1,n) = eff(n);
    T1map(n,1) = (1/r1(n))*1000;
end
disp('Fitting T1...Finished.')

disp('Fitting T2...')
TEdata = reshape(TEdata,[size(TEdata,1)*size(TEdata,2), size(TEdata,3)]);
fittedT2map = zeros(1,size(TEdata,1));
T2model = @(a,T2,c, x) a*exp(-x*(1/T2)) + c;
for n = 1:size(TEdata,1)
    [fittedT2Curve] = fit(TE',TEdata(n,:)',T2model,'Upper',[5000 5000 2500],'Lower',[0 0 0],'StartPoint',[1000 100 100]);
    T2map(n,1) = fittedT2Curve.T2;
end
disp('Fitting T2...Finished.')

fittedT2map = reshape(T2map,64,64);
fittedT1map = reshape(T1map,64,64);

%% 4. Plot and Save Results
figure,
subplot 121
imagesc(fittedT1map)
c = colorbar;
caxis([0 1000])

subplot 122
imagesc(fittedT2map)
c = colorbar;
caxis([0 1000])


%%
matlab2tikz([savingdir,'/MAT-files/customPhantom2goldstandardfittedT1T2',date])

