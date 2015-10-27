function plotSimulatedSignal(M, Mxy, imageTimes)

figure
hold on
plot(M(1,:))
plot(M(2,:))
plot(M(3,:))
plot(abs(complex(M(1,:),M(2,:))),'--k')
plot(imageTimes, Mxy,'+')
xlabel 'Time (ms)'
ylabel 'Magnetisation'
legend 'Mx' 'My' 'Mz' 'Mtransverse'
set(gca,'fontsize',20)

end