
figure; 
plot(SNR,tSNR,'rx'); axis equal; xlim([0 max(SNR)*1.1]);
ylim([0 max(SNR)*1.1]); hold on;
plot([0 1]*max(SNR),[0 1]*max(SNR),'k--'); grid on;
xlabel 'SNR'; ylabel 'tSNR';
legend('Measured values','Equality')