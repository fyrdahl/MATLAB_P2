function plot_simulated_signal(M, Mxy, imageTimes)

figure('name',['Simulated Magnetisation. ',datestr(now,30)])
hold on
plot(M(1,:))
plot(M(2,:))
plot(M(3,:))
plot(abs(complex(M(1,:),M(2,:))),'--k')
for i = 1:size(Mxy,1)
plot(imageTimes, Mxy(i,:),'+')
end
xlabel 'Time [ms]'
ylabel 'Magnetisation [a.u.]'
legend('Mx', 'My', 'Mz' ,'Mtransverse','Image Signal','location','best')
set(gca,'fontsize',20)

end