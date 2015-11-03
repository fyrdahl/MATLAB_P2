function plotSimulatedSignal(M, Mxy, imageTimes,offsetListNum)

figure('name',['Simulated Magnetisation. Offset List:',num2str(offsetListNum)])
hold on
plot(M(1,:))
plot(M(2,:))
plot(M(3,:))
plot(abs(complex(M(1,:),M(2,:))),'--k')
plot(imageTimes, Mxy,'+')
xlabel 'Time [ms]'
ylabel 'Magnetisation [a.u.]'
legend('Mx', 'My', 'Mz' ,'Mtransverse','Image Signal','location','best')
set(gca,'fontsize',20)

end