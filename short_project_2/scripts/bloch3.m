%Using the stanford tutorial
T = 1000;
dt = 1;
TE = 50;
TR = 1000;
N = ceil(T/dt) + 1;
M = zeros(3,N);
M(:,1) = [0, 0, 1]';
M0 = 1; %fully relaxed equilibrium
T1 = 600;
T2 = 100; %ms

N1 = round(TE/2/dt);
N2 = round((TR-TE/2)/dt);


el = 1;

df = [-50, -40,-30,-20,-10, 0,10,20, 30, 40, 50];
clear Msig
for f = 1:length(df)
[A, B] =  freeprecess(dt, T1, T2, df(f));

M(:,el+1) = Rot_y(90)*M(:,1);
for i = 3:(N1+1)
    M(:,i) = A*M(:,i-1) + B;
end
M(:,N1+2) = Rot_x(180)*M(:,N1+1) + B;
for i = 2:N2-1
    M(:,i+N1+1) = A*M(:,i+N1) + B;
end
Msig(:,f) = M(1,:)+i*M(2,:);
end
figure
time = [0:N-1]*dt;
hold on
%plot(time,M(1,:),'b-',time,M(2,:),'r--',time,M(3,:),'g-.');
plot(time,abs(mean(Msig,2)),'b-');
%plot(time,abs(Msig(:,1)),'k.');
legend('M_x','M_y','M_z');
xlabel('Time (ms)');
ylabel('Magnetization');
%axis([min(time) max(time) -1 1]);
grid on;


