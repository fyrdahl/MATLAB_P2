
v = VideoWriter('dictionaryMov.avi');

v.FrameRate = 2;

open(v)


figure;

d = squeeze(FPdata(25,25,:));
d =    d(:)/d(1);

for n = 1:60:size(dictionary,1)
    
    title(['Dictionary Entry:',num2str(n),' of ',num2str(size(dictionary,1))])
    dict = dictionary(n,:)/dictionary(n,1);
    plot(dict,'b-')
   % set(gca,'FontSize',24)
    ylabel(['Normalised Signal Intensity'])
    xlabel(['Image Index'])
    ylim([ 0 2.2])
    hold on
    drawnow
    frame = getframe;
    writeVideo(v,frame);
    
    
end
close(v);

plot(d,'r-o','lineWidth',2)
set(gca,'FontSize',24)
ylabel(['Normalised Intensity'])
xlabel(['Image Index'])

legend('Dictionary Signals','Data Example','Location','Best')
