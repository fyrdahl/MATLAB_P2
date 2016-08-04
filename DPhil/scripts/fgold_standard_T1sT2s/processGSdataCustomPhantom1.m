%% Gold Standard T1 and T2 maps of Custom Phantom 1 

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
ID = 'MR_PHYSICS_359';
date = '20150819';
phantom = 'Jack1';
TIprotocol = 'IR_ep2d_se_T1_';
TEprotocol = 'pj_ep2d_se_TE';
TI =  [90:30:150, 200:50:400, 500:100:1000, 1500, 2000:1000:5000];
TE = [32,40:10:90,100:20:160];
TIseriesStart = 4;
TEseriesStart = 23;

for n = 1:6
TIimage(:,:,n) = nifti2mat([savingdir,'/raw_data/NIfTI/',date,'_',ID,'/TIgroup/',date,'_00',num2str(TIseriesStart+(n-1)),'_',TIprotocol,num2str(TI(n)),'ms.nii']);
end
for n = 7:numel(TI)
TIimage(:,:,n) = nifti2mat([savingdir,'/raw_data/NIfTI/',date,'_',ID,'/TIgroup/',date,'_0',num2str(TIseriesStart+(n-1)),'_',TIprotocol,num2str(TI(n)),'ms.nii']);
end

for n = 1:numel(TE)
TEimage(:,:,n) = nifti2mat([savingdir,'/raw_data/NIfTI/',date,'_',ID,'/TEgroup/',date,'_0',num2str(TEseriesStart+(n-1)),'_',TEprotocol,'_',num2str(TE(n)),'ms.nii']);
end

%%

offsetList(:,1) = 2000;
offsetList(:,2) = [32,40:10:90,100:20:160];
offsetList(:,3) = 90;
offsetList(:,3) = 180;

plot_sim_comparison(TEimage,coords, 10, T2, offsetList,nRepeats,df)


%% Fit T1 and T2 curves to Gold Standard images
TIdata = reshape(TIimage,[size(TIimage,1)*size(TIimage,2),size(TIimage,3)]);
opts = struct('debug',0,'fiteff',1);

disp('Started T1 fits...')
for n = 1:size(TIdata,1)
    n
    [pd(n) r1(n) eff(n)] = qmap_t1_fit_ir(TIdata(n,:), 0.001*TI,opts);
    BetaMap(1,n) = eff(n);
    T1map(n,1) = (1/r1(n))*1000;
end
disp('Finished T1 fits.')

disp('Started T2 fits...')
TEdata = reshape(TEimage,[size(TEimage,1)*size(TEimage,2), size(TEimage,3)]);
fittedT2map = zeros(1,size(TEdata,1));
T2model = @(a,T2,c, x) a*exp(-x*(1/T2)) + c;
for n = 1:size(TEdata,1)
    n
    [fittedT2Curve] = fit(TE',TEdata(n,:)',T2model,'Upper',[5000 5000 2500],'Lower',[0 0 0],'StartPoint',[1000 100 100]);
    T2map(n,1) = fittedT2Curve.T2;
end
disp('Finsihed T2 fits.')

fittedT2map = reshape(T2map,64,64);
fittedT1map = reshape(T1map,64,64);

%%
figure, 
subplot 121
imagesc(fittedT1map)
c = colorbar;
caxis([0 3500])

subplot 122
imagesc(fittedT2map)
c = colorbar;
caxis([0 3000])

