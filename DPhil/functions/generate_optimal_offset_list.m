% generate optimal offset list
% calculate signal difference for a range of T1 and T2, using a range of
% different TEs and TRs and flip angles (FA1 and FA2).
t=1:2000;
T1 = 282.3;
T2 = 214.8;

signal = exp(-t/T2);
max(t(signal >= 0.95*signal(1)))
min(t(signal <= 0.05*signal(1)))
n = unifrnd(max(t(signal >= 0.95*signal(1))),min(t(signal <= 0.05*signal(1))),24,1);

plot(signal)
hold on
plot(round(n),signal(round(n)),'*')

Mz = 1 - exp(-t/T1);
min(t(Mz >= 0.95*1))
max(t(Mz <= 0.05*1))

offsetList = generate_offset_list([max(t(Mz <= 0.05*1)) 12 90 180],[min(t(Mz >= 0.95*1)) 645 90 180],24,1);
