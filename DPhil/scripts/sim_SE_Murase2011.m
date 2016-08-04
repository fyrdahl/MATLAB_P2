%function [absMsig, M] = sim_SE_Murase2011(T1, T2,offsetList, df)

T1 = 3000;
T2 = 300;

offsets(1:5,1) = 2000;
offsets(1:5,2) = 1000;
offsets(1:5,3) = 90;
offsets(1:5,4) = 180;
[Msig, M] = sim_SE_bernstein(T1*1000, T2*1000, 0, 0, offsets,1,0);



%M = [Mx My Mz 1]'
T1 = 280/1000; %seconds
T2 = 214/1000; %seconds
Mz0 = 1;
M = [0 0 Mz0 1]';
B1 = 3; %RF power
R1 = 1/T1;
R2 = 1/T2;
omega1 = 42580000/B1;%42.58 MHz/Tesla, larmor frequency
dOmega = 0; %Hz

    A = [ -R2 dOmega 0 0; ...
    -dOmega -R2 0 0; ...
    0 -0 -R1 R1*Mz0; ...
    0 0 0 0];

for n = 1:5

    % FA1
M = rot_y(offsets(n,3))*M(1:3);

    % to TE/2
M = expm(A*(offsets(n,2)/2/1000))*[M;1];

    % to FA2
M = rot_y(offsets(n,4))*M(1:3);

    % to TE
M = expm(A*((offsets(n,2)/2)/1000))*[M;1]

    % to TR
M = expm(A*((offsets(n,1)-offsets(1,2))/1000))*M;


end



