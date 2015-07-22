addpath(genpath('/Applications/fsl/'))
addpath(genpath('/usr/local/fsl'))
addpath(genpath('/usr/local/fsl/bin'))
addpath(genpath('/Users/jallen/Documents/MATLAB'))
%% Read images
rawTEImages = dir('/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/20150714_100527ep2dseTE*');
rawTIImages = dir('/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/20150714_100527IRep2dseTI*');

TE = [32:10:72,92:20:152,192:50:292,372,400];
TI = [35,85,135,185,235,285,335,385,435,485,585,685,785,885,985];

TEimages = zeros(64,64,TE(end));
TIimages = zeros(64,64,TI(end));


for TE_ind = 1:size(rawTEImages,1)  
    for i = 1:size(rawTEImages,1)
        i
        if strfind(rawTEImages(i).name(1,:),['TE',num2str(TE(TE_ind))]) ~= 0
            filename = ['/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/', rawTEImages(i).name(1,:)];
            TEimages(:,:,TE(TE_ind)) = read_avw(filename);
            
        end
    end
end

for TI_ind = 1:size(rawTIImages,1)   
    for i = 1:size(rawTIImages,1)
        if strfind(rawTIImages(i).name(1,:),['TI',num2str(TI(TI_ind))]) ~= 0
            filename = ['/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/', rawTIImages(i).name(1,:)];
            TIimages(:,:,TI(TI_ind)) = read_avw(filename);        
        end
    end
end
%% Calculate T1 and T2
% set the position and size of the ROI from which to calculate the mean
ROI_dim = [20:25; 20:25]

% figure
% imagesc(TEimages(:,:,1))
% hold on
% rectangle('Position',[ROI_dim(1,1),ROI_dim(2,1),size(ROI_dim,2),size(ROI_dim,2)])

for i = TE
ROI = (TEimages(ROI_dim(1,:),ROI_dim(2,:),i))
TEmeans(i) = (mean(ROI(:)))
end
for i = TI
ROI = (TIimages(ROI_dim(1,:),ROI_dim(2,:),i))
TImeans(i) = (mean(ROI(:)))
end

figure('name','Varying TE')
plot(TE,TEmeans(1,TE),'*')
ylim([0 3500]);
hold on
exp_fit1 = fit([TE]',(TEmeans(1,TE))','exp1');
plot(TE,exp_fit1.a*exp((exp_fit1.b)*TE))
xlabel 'TE'
ylabel 'Signal'
T2 = - 1/exp_fit1.b

newTImeans = TImeans;
newTImeans(1,TI(1:4)) = -newTImeans(1,TI(1:4));
newTImeans = -newTImeans;
newTImeans = newTImeans + 3500;
figure('name','Varying TI')
hold on
plot(TI, TImeans(1,TI),'*')
plot(TI, newTImeans(1,TI),'*')
plot([0 1000],[0 0])
exp_fit2 =fit(TI',(newTImeans(1,TI))','exp1');
% plot(TI,(exp_fit2.a)*exp((exp_fit2.b)*TI));
plot(TI,(exp_fit2.a-3500)*(1-2*exp((exp_fit2.b)*TI)));
xlabel 'TI'
ylabel 'Signal'
T1 = -1/exp_fit2.b % ms
% legend ([num2str(M),'(1-2*exp(,',num2str(A),'*TI'])
%%
run('blochSim')

