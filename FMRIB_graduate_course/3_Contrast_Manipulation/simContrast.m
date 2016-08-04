% function simContrast
% clear; ca

% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

if ~exist('TR','var') || ~exist('flipAngle','var') || ~exist('TE','var')
    disp('ERROR: TR, flipAngle and TE must all be defined')
    return
end

% 1- WM
% 2- GM
% 3- CSF
if ~exist('fieldStrength','var'),    fieldStrength=1.5; end
if fieldStrength==3
    T1vals = [830 1330 4000];
    T2vals = [80 110 2000];
else
    fieldStrength=1.5;
    T1vals = [600 900 3500]; % T1 in ms
    T2vals = [80 100 2000];  % T2 in ms    
end

tColours = {'b',[0 .7 0],'r'};

load segBrain

% TR = 5000e-3;
% TE = 90e-3;
if ~exist('useInversion','var'), useInversion = 0; end
% useInversion = 1;
if (useInversion) && ~exist('TI','var')
    disp(['ERROR: if inversion is chosen, TI must be set']);
    return
end
% TI = 600;
flip = flipAngle;

maxTE = min(TR,200);

if TE>maxTE
    disp(['ERROR: TE may not be greater than TR or 200ms (whichever is smaller)'])
    return
end

if useInversion && (TI>=TR)
    disp(['ERROR: TI may not be greater than (or equal to) TR'])
    return
end

figure
% set(gcf,'Position',[    50   164   592   945])
set(gcf,'Position',[    50   512   937   597])
% axHeight = .22; topY = .73; midY = .4; botY = .05; thinWidth = .08;
% axLeft = .12; mainWidth = .6876; loff = .1;
axHeight = .35; topY = .58; midY = .1; thinWidth = .05; axLeft = .09; mainWidth=.5;loff=.02;
subT1 = subplot('Position',[    axLeft    topY    mainWidth    axHeight]);
hold on; grid on; xlabel('Time since last excitation pulse (ms)')
if useInversion
    title(['{\bf T1-weighting:} TR = ' num2str(TR) 'ms, Flip Angle = ' num2str(flip) ' degrees, TI = ' num2str(TI) 'ms'])
else
    title(['{\bf T1-weighting:} TR = ' num2str(TR) 'ms, Flip Angle = ' num2str(flip) ' degrees'])
end
ylabel('Relative M_z')
    
subT1w = axes('Position',[    mainWidth+axLeft+loff   topY    thinWidth    axHeight]);
set(gca,'xtick',[],'ytick',[],'box','off')
hold on; axis([0 2 0 1]);
title('T1-weighting')
subT2 = axes('Position',[   axLeft    midY   mainWidth    axHeight]);
hold on; grid on; title(['{\bf T2-weighting:} TE = ' num2str(TE) 'ms'])
xlabel('Time since excitation pulse (ms)')
ylabel('Relative M_{xy}')
axis([0 maxTE 0 1.05])
subT2w = subplot('Position',[    mainWidth+axLeft+loff    midY    thinWidth    axHeight]);
hold on; axis([0 2 0 1]);
set(gca,'xtick',[],'ytick',[],'box','off')
title('T2-weighting')

subText = axes('Position',[    0.0164    0.9397    0.1010    0.0457]);
axis off
axis([0 1 0 1])
text(0,1,['{\bfField Strength:}' num2str(fieldStrength) 'T'])

% 
% subT1w2 = subplot('Position',[    axLeft-.02    botY   thinWidth
% axHeight]);
% set(gca,'xtick',[],'ytick',[],'box','off')
% ylabel('Relative Signal')
% hold on; axis([0 2 0 1]);title('\bf T1')
% subT2w2 = subplot('Position',[    0.2-.02   botY    thinWidth   axHeight]);
% hold on; axis([0 2 0 1]); title('\bf T2')
% set(gca,'xtick',[],'ytick',[],'box','off')

% subT1T2w = subplot('Position',[    0.34   botY   .18    axHeight]);
subT1T2w = subplot('Position',[ 0.7278    0.6913    0.1987    0.2381]);

hold on; axis([0 2 0 1.1]); title({'\bf Combined contrast','(normalised)'})
set(gca,'xtick',[.5 1 1.5],'ytick',[],'xticklabel',{'WM','GM','CSF'})
ylabel('Relative Signal')


% subIm = subplot('Position',[ .55 botY .4 axHeight]);
subIm = subplot('Position',[    0.6931    0.0688    0.2517    0.4948]);


lw = 10;

if useInversion
    M_minus = (1 - exp(-(TR-TI)./T1vals) + (1 - exp(-TI./T1vals)).*cos(flip*pi/180).*exp(-(TR-TI)./T1vals))./(1 + exp(-TR./T1vals)*cos(flip*pi/180));
    MTI_minus = ((1 - exp(-TI./T1vals)) - M_minus.*exp(-TI./T1vals));
    sigT1 = abs(MTI_minus*sin(flip*pi/180));
    t = linspace(0,TR,500);
    iTI = round(interp1(t,1:length(t),TI));
    M = zeros(length(t),length(T1vals));
    for iT1 = 1:length(T1vals)
        M(1:iTI,iT1) = (1 - exp(-t(1:iTI)/T1vals(iT1))) - M_minus(iT1)*exp(-t(1:iTI)/T1vals(iT1));
        M(iTI+1:end,iT1) = (1 - exp(-(t(iTI+1:end)-TI)/T1vals(iT1))) + MTI_minus(iT1)*cos(flip*pi/180)*exp(-(t(iTI+1:end)-TI)/T1vals(iT1));
    end
 
else
    M_minus = (1-exp(-TR./T1vals))./(1-cos(flip*pi/180)*exp(-TR./T1vals));
    sigT1 = M_minus*sin(flip*pi/180);
    t = linspace(0,TR,500);
    M = zeros(length(t),length(T1vals));
    for iT1 = 1:length(T1vals)
        M(:,iT1) = 1 - exp(-t/T1vals(iT1)) + M_minus(iT1)*cos(flip*pi/180)*exp(-t/T1vals(iT1));
    end


end
tt = [-t(125:-1:1) t  t(1:125)+TR];
for iLine = 1:length(T1vals)
    subplot(subT1)
    plot(tt,[M(end-124:end,iLine); M(:,iLine); M(1:125,iLine)],'linewidth',2,'color',tColours{iLine});
    axis([tt(1) tt(end) min( min(M(:)) ,0) 1.05*max(M(:))])
        
    subplot(subT1w)
    line([0 0]+iLine*.5,[0 sigT1(iLine)],'linewidth',lw,'color',tColours{iLine})
%     subplot(subT1w2)
%     line([0 0]+iLine*.5,[0 sigT1(iLine)],'linewidth',lw,'color',tColours{iLine})
end
subplot(subT1)

if ~useInversion
    line([0 0], [-1 1], 'linewidth',2,'color','k','linestyle','--')
    line([TR TR], [-1 1], 'linewidth',2,'color','k','linestyle','--')
else
    line([0 0], [-1 1], 'linewidth',4,'color','k','linestyle',':')
    line([TR TR], [-1 1], 'linewidth',4,'color','k','linestyle',':')
    line([TI TI], [-1 1], 'linewidth',2,'color','k','linestyle',':')
    xlabel('Time since last inversion pulse (ms)')
end


t = linspace(0,maxTE,100);
sigT2 = exp(-TE./T2vals);
sigT1T2 = sigT1.*sigT2;
subplot(subT2)
line([TE TE],[0 1],'linewidth',2,'color','k','linestyle','--')
for iLine = 1:length(T2vals)
    subplot(subT2)
    plot(t,exp(-t/T2vals(iLine)),'linewidth',2,'color',tColours{iLine}); 
    
    subplot(subT2w)
    line([0 0]+iLine*.5,[0 sigT2(iLine)],'linewidth',lw,'color',tColours{iLine})
   
%     subplot(subT2w2)
%     line([0 0]+iLine*.5,[0 sigT2(iLine)],'linewidth',lw,'color',tColours{iLine})
    
    subplot(subT1T2w)
    line([0 0]+iLine*.5,[0 sigT1T2(iLine)/max(sigT1T2)],'linewidth',20,'color',tColours{iLine})
end

fontScale(1.3)
bIm = zeros(size(segBrain));
for iLine = 1:length(T1vals)
    bIm(segBrain==iLine) = sigT1T2(iLine)/max(sigT1T2);
end
subplot(subIm)
noiseIm = randn(size(segBrain))*.005/max(sigT1T2);
imagesc(abs(bIm+noiseIm),[0 1]); axis equal tight off; colormap(gray)
% title({['rel. GM/WM contrast: ' abs(num2str(abs(sigT1T2(2)-sigT1T2(1))/sigT1T2(1),3))],...
title({    ['{\bf Total acq. time:} ' num2str(round(192*128*TR/1000/6)/10) ' mins'],...
    ['{\bf SNR:} ' num2str(max(sigT1T2)*200,3) ', {\bf GM/WM CNR}: ' num2str(abs((sigT1T2(2)-sigT1T2(1))*200),3)] })

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')

