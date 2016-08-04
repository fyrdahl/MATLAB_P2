
T1 = 500;
T2 = 400;
omega = 42.5781*(2*pi); % MHz/T
domega = 0; %MHz/T

FA = 90;
Mzeq = 1;
M0 = [ 0 0 Mzeq 1]'


%rotation matrix 1
R1(1,1:4) = ([1, 0 , 0, 0]);
R1(2,1:4) = ([0, cosd(FA), -sind(FA), 0]);
R1(3,1:4) = ([0, sind(FA) , cosd(FA), 0]);
R1(3,1:4) = ([0, sind(FA) , cosd(FA), 0]);
R1(4,:) = ([0, 0 , 0, 1]);

n = 1;
A = [ -1/T2, domega, 0, 0;...
    -domega, -1/T2, omega, 0;...
    0,  -omega, -1/T1, (1/T1)*Mzeq;...
    0,  0,  0,  0];
R1*M0;
expm(A*tau);
M = expm(A*tau)*R1*M0;
n = n + 1;

%%

df = 0;
M = [0 0 1]';

T = 100;

phi = 2*pi*df*T/1000;% Resonant precession, radians.
rot_z(phi);
%evolution matrix
E = [exp(-T/T1) 0 0;
     0  exp(-T/T2) 0;
     0          0 1-exp(-T/T1)];








% Flip 2, a rotation about the x axis
M(:,n) = rot_x(parameters(n,4))*M(:,n);

% Magnetisation evolution over time until TE
[Afp]=freeprecess(((parameters(n,2)+TEmin)/2),T1,T2,df);
M(1,n) = M(1,n)*Afp(1,1);
M(2,n) = M(2,n)*Afp(2,2);
M(3,n) = Mzeq - (Mzeq - M(3,n))*Afp(3,3);
% Collect Signal
Msig(n) = complex(M(1,n),M(2,n));
absMsig(n) = abs(Msig(n));

% if the phase of the RF pulse changes, we need Fnx Fstarnx Fny Fstarny
%EchoAmplitudes(n) = sqrt(Fstarnx^2 + Fstarny^2)

% At the end of the TR, set the current Mz to be the initial Mz for the next TR

% TRmin + nSlices*(TEmin + TEadd) + TRadd
[Afp]=freeprecess(parameters(n,2) + TRmin + parameters(n,1) - (parameters(n,2)+TEmin),T1,T2,df);
M(1,n) = M(1,n)*Afp(1,1);
M(2,n) = M(2,n)*Afp(2,2);
M(3,n) = Mzeq - (Mzeq - M(3,n))*Afp(3,3);
M0 = [ 0;  0; M(3,n)]; % ignore the remaining transverse signal


%% plotting the signal evolution
if strcmp(plotFlag,'plot')
    figure
    hold on
    plot(M(1,:))
    plot(M(2,:))
    plot(M(3,:))
    plot(TRindex, Msig,'+')
    xlabel 'Time (ms)'
    ylabel 'Magnetisation'
    legend 'Mx' 'My' 'Mz' 'Mtransverse'
    set(gca,'fontsize',20)
end


