
t = [1:12000];

T1(1) = 282;
T1(2) = 1300;
T1(3) = 700;
T1(4) = 3000;

T2(1) = 214;
T2(2) = 100;
T2(3) = 60;
T2(4) = 2500;

TE = [25:10:175, 200:100:1000];
TI = [100:50:600,700:100:1000, 1250:250:4500];

for n = 1:numel(T1)
Mxy(n,:) = exp(-t/T2(n));
Mz(n,:) = 1-exp(-t/T1(n));
end

figure,
title 'Homogeneous Phantom'
hold on
for n = 1
    plot(Mxy(n,:))
    plot(Mz(n,:))
    plot(TE,exp(-TE/T2(n)),'*')
    plot(TI,(1-exp(-TI/T1(n))),'+')
    
end

xlabel 'Time [ms]'
ylabel 'Magnetisation'

figure,
title 'Custom Phantom 2'
hold on
for n = 2:numel(T2)
    plot(Mxy(n,:))
    plot(Mz(n,:))
    plot(TE,exp(-TE/T2(n)),'*')
    plot(TI,(1-exp(-TI/T1(n))),'+')
    
end

xlabel 'Time [ms]'
ylabel 'Magnetisation'