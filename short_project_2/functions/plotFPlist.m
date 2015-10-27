function plotFPlist(fingerprintLists)
list = squeeze(fingerprintLists);

figure
subplot 221
plot(list(:,1),'-*')
xlabel 'Index'
ylabel 'Repetition Time (TR) [ms]'
subplot 222
plot(list(:,2),'-*')
xlabel 'Index'
ylabel 'Echo Time (TE) [ms]'
subplot 223
plot(list(:,3),'-*')
xlabel 'Index'
ylabel 'Flip Angle 1 [degrees]'
subplot 224
plot(list(:,4),'-*')
xlabel 'Index'
ylabel 'Flip Angle 2 [degrees]'
end