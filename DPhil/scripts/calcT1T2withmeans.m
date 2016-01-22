%% Calculate T1 and T2: with means

% set the position and size of the ROI from which to calculate the mean
ROI_dim = [20:25; 20:25];
ROI_dim = [23:26; 27:30];

for i = TE
    ROI = (TEimages(ROI_dim(1,:),ROI_dim(2,:),i));
    TEmeans(i) = (mean(ROI(:)));
end
for i = TI
    ROI = (TIimages(ROI_dim(1,:),ROI_dim(2,:),i));
    TImeans(i) = (mean(ROI(:)));
end

figure('name','Varying TE')
plot(TE,TEmeans(1,TE),'*')
ylim([0 3500]);
hold on
exp_fit1 = fit([TE(2:end)]',(TEmeans(1,TE(2:end)))','exp1');
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

plot(TI,(exp_fit2.a-3500)*(1-2*exp((exp_fit2.b)*TI)));
xlabel 'TI'
ylabel 'Signal'
T1 = -1/exp_fit2.b ; % ms