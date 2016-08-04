% Single calculation with EPG calculations

% Magnetisation vector, including dephased and rephased populations
T1 = 500; %ms
T2 = 300; %ms

ETL = 2
MPSV_0 = zeros(3*ETL,1);
MPSV_0(1) = 1;


FA = 180;

% Rotation matrix to apply to all states
R_0 = [cosd(FA/2).^2  ,   sind(FA/2).^2 ,  -1i*sind(FA);...
      sind(FA/2).^2 ,    cosd(FA/2).^2,  1i*sind(FA) ;...
     -0.5*1i*sind(FA),  0.5*1i*sind(FA),    cosd(FA)];
R = blkdiag(R_0,R_0)
for n = 3:ETL
R = blkdiag(R,R_0)
end

% Transition Matrix, for between RF pulse
T_0 = rot_x(FA);
T_1 = [cosd(FA/2).^2  ,   sind(FA/2).^2 ,  sind(FA) ;...
      sind(FA/2).^2 ,    cosd(FA/2).^2,    0 ;...
     0.5*sind(FA),   -0.5*sind(FA) ,    cosd(FA)]
T = blkdiag(T_0,T_1)





%Prepare relation matrix E
TE = 100; %Echo Time [ms]
E = blkdiag(exp(-TE/(2*T1)),exp(-TE/(2*T2)));
E = blkdiag(E,E)

M(n) = MPSV(1)
%% Signal evolution

    % 1. RF pulse 1, a rotation about x-axis
M = rot_x(90)*M;









    % incorporate EPG algorithm (see Hennig '91)
    F = complex(M(1,n),M(2,n));
    Fstar = complex(M(1,n),-M(2,n));
    
    % Magnetisation evolution over time until TE/2
    [Afp]=freeprecess(((parameters(n,2)+TEmin)/2),Trefoc,T2,df);
    M(1,n) = M(1,n)*Afp(1,1);
    M(2,n) = M(2,n)*Afp(2,2);
    M(3,n) = Mzeq - (Mzeq - M(3,n))*Afp(3,3);
    
    % Flip 2, a rotation about the x axis
    M(:,n) = rot_x(parameters(n,4))*M(:,n);
    
    % Magnetisation evolution over time until TE
    [Afp]=freeprecess(((parameters(n,2)+TEmin)/2),Trefoc,T2,df);
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
    [Afp]=freeprecess(parameters(n,2) + TRmin + parameters(n,1) - (parameters(n,2)+TEmin),Trefoc,T2,df);
    M(1,n) = M(1,n)*Afp(1,1);
    M(2,n) = M(2,n)*Afp(2,2);
    M(3,n) = Mzeq - (Mzeq - M(3,n))*Afp(3,3);   
    M0 = [ 0;  0; M(3,n)]; % ignore the remaining transverse signal
end

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

end
