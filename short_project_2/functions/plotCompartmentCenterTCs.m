function plotCompartmentCenterTCs(compartmentCenters,TEimages, TIimages, TE, TI)


    
figure('name', 'Inversion Recovery')
for i = 1:6
    hold on
    plot(TE(2:end),squeeze(TEimages(compartmentCenters(i,1,1),compartmentCenters(i,2,1),TE(2:end))),'-*')
end
ylabel 'Signal [a.u.]'
xlabel 'Echo Time (TE) [ms]'

figure('name', 'Spin Echo')
for i = 1:6
    hold on
    plot(TI(2:end), squeeze(TIimages(compartmentCenters(i,1,2),compartmentCenters(i,2,2),TI(2:end))),'-*')
  
end
ylabel 'Signal [a.u.]'
xlabel 'Inversion Time (TI) [ms]'


end