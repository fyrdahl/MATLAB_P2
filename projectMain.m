addpath(genpath('/Applications/fsl/'))
addpath(genpath('/usr/local/fsl'))
addpath(genpath('/usr/local/fsl/bin'))
addpath(genpath('/Users/jallen/Documents/MATLAB'))
%% Read images
rawTEImages = dir('/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/20150714_100527ep2dseTE*');
rawTIImages = dir('/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/20150714_100527IRep2dseTI*');

rawTEImages = dir('/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/20150806_125113ep2dseTE*');
rawTIImages = dir('/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/20150806_125113IRep2dseTI*');


%%
TE = [32:10:72,92:20:152,192:50:292,372,400];
TI = [35,85,135,185,235,285,335,385,435,485,585,685,785,885,985];

TE = [32:10:72,92:20:152,192:50:292,372,400];
TI = [35,85,135,185,235,285,335,385,435,485,585,685,785,885,985];



for i = 1:numel(rawTEImages)
   
    name = rawTEImages(i).name
    if name(1,27)==['s']
        tmpTE(i,1:2) = name(1,24:25)
    else 
        tmpTE(i,1:3) = name(1,24:26)
    end 
   TE(i) = str2num(tmpTE(i,:))

end
TE = unique(sort(TE))

for i = 1:33
    name = rawTIImages(i).name
    if name(1,28)==['m']
        tmpTI(i,1:2) = name(1,26:27)
    end
    if name(1,29)==['m']
        tmpTI(i,1:3) = name(1,26:28)
    end
    if name(1,30)==['m']
        tmpTI(i,1:4) = name(1,26:29)
    end 
  
   TI(i) = str2num(tmpTI(i,:))

end
TI = unique(sort(TI))

%%
TEimages = zeros(64,64,TE(end));
TIimages = zeros(64,64,TI(end));

for TE_ind = 1:numel(TE)  
  
    for i = [1:13,17:size(rawTEImages,1)]
       
        if strfind(rawTEImages(i).name(1,:),['TE',num2str(TE(TE_ind))]) > 0
            filename = ['/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/', rawTEImages(i).name(1,:)]
            TEimages(:,:,TE(TE_ind)) = read_avw(filename);
            test = TE(TE_ind)
            pause
        end
    end
    
end

figure
for i = 1:numel(TE)
    colormap gray
imagesc(TEimages(:,:,TE(i)))

pause
end

for TI_ind = 1:size(rawTIImages,1) 
    for i = [1:17,20:size(rawTIImages,1)]
        
        if strfind(rawTIImages(i).name(1,:),['TI',num2str(TI(TI_ind))]) > 0
            filename = ['/Users/jallen/Documents/MATLAB/short_project_2/nifti_data/', rawTIImages(i).name(1,:)];
      TIimages(:,:,TI(TI_ind)) = read_avw(filename);
        end
    end   
end


%% Calculate T1 and T2
% set the position and size of the ROI from which to calculate the mean
ROI_dim = [20:25; 20:25];

% figure
% imagesc(TEimages(:,:,1))
% hold on
% rectangle('Position',[ROI_dim(1,1),ROI_dim(2,1),size(ROI_dim,2),size(ROI_dim,2)])

for i = TE
ROI = (TEimages(ROI_dim(1,:),ROI_dim(2,:),i));
TEmeans(i) = (mean(ROI(:)));
end
for i = TI
ROI = (TIimages(ROI_dim(1,:),ROI_dim(2,:),i));
TImeans(i) = (mean(ROI(:)));
end

 fitted = zeros(64,64)

for r = 1:64
    r
for c = 1:64
    c;
   fittedModel = fit(TE',squeeze(TEimages(r,c,TE)),'exp1');
   fitted(r,c) = fittedModel.b;
end
end


figure('name','Varying TE')
plot(TE,TEmeans(1,TE),'*')
ylim([0 3500]);
hold on
exp_fit1 = fit([TE]',(TEmeans(1,TE))','exp1');
plot(TE,exp_fit1.a*exp((exp_fit1.b)*TE));
xlabel 'TE'
ylabel 'Signal'
T2 = - 1/exp_fit1.b

newTImeans = TImeans;
newTImeans(1,TI(1:4)) = -newTImeans(1,TI(1:4));
newTImeans = -newTImeans;
newTImeans = newTImeans + 3500;
figure('name','Varying TI')
hold on
plot(TI, TImeans(1,TI),'*');
plot(TI, newTImeans(1,TI),'*');
plot([0 1000],[0 0]);
exp_fit2 =fit(TI',(newTImeans(1,TI))','exp1');
% plot(TI,(exp_fit2.a)*exp((exp_fit2.b)*TI));
plot(TI,(exp_fit2.a-3500)*(1-2*exp((exp_fit2.b)*TI)));
xlabel 'TI'
ylabel 'Signal'
T1 = -1/exp_fit2.b ; % ms
% legend ([num2str(M),'(1-2*exp(,',num2str(A),'*TI'])
%%
run('blochSim')

%%
t0 = 0
for i = 1:24
    
   
    
A(i) = TEmeans(TE(i))/TEmeans(TE(1))
C(i) = imag(Mtransverse(t + TE(i) ))/imag(Mtransverse(1));

t = t + newTR(TE(i))
end

for i = TI
B(i) = TImeans(i)/TImeans(TI(end))
end

for i = 1:numel(Mtransverse)
C(i) = imag(Mtransverse(i))/imag(Mtransverse(1));
end
