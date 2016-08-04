function M = show_traj2(Gx,Gy)

% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

  fig=figure;    
  set(fig,'Position',[70 160 800 720]);
  fprintf('Hit enter to continue...\n');
  pause
  for j=1:length(Gx),
    plot_index(Gx,Gy,j);
    M(j)=getframe(fig);
  end;


function [] = plot_index(Gx,Gy,j)

    k=cumsum(Gx+i*(Gy+1e-6));
    m=max([abs(Gx) abs(Gy)]);
    kmax=1.05*max(max(real(k)),max(imag(k)));
    t=1:length(k);
  
    subplot('position',[0.03 0.57 0.45 0.35]);
    hold off; plot(Gx,'LineWidth',2); hold on; 
    plot([t(j) t(j)],[-m m],'r-','LineWidth',2); 
    plot([t(1) t(end)],[0 0],'k-','LineWidth',2);
    axis off; tt=text(0,1.1*m,'Gx'); set(tt,'FontSize',18);
    
    subplot('position',[0.03 0.12 0.45 0.35]);
    hold off; plot(Gy,'LineWidth',2); hold on; 
    plot([t(j) t(j)],[-m m],'r-','LineWidth',2); 
    plot([t(1) t(end)],[0 0],'k-','LineWidth',2);
    axis off; tt=text(0,1.1*m,'Gy'); set(tt,'FontSize',18);

    subplot('position',[0.5 0.25 0.45 0.55]);
    hold off; plot(k(1:j),'LineWidth',2); hold on;
    plot([-kmax kmax],[0 0],'k-','LineWidth',2); 
    plot([0 0],[-kmax kmax],'k-','LineWidth',2);
    axis([-kmax kmax -kmax kmax]); 
    plot(k(j),'rx','LineWidth',2);
    axis off; tt=text(-kmax,1.1*kmax,'k-space'); set(tt,'FontSize',18);
    tt=text(-1.1*kmax, 0, 'kx');  set(tt,'FontSize',18);
    tt=text(0,-1.1*kmax,'ky'); set(tt,'FontSize',18);

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')
