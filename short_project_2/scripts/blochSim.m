%% Bloch simulation
%JALLEN
%%

clear M
clear Mtransverse
T1 = 292; % ms
T2 = 214; % ms
tau = [500; 500; 500; 500; 500;];
TE = 2*tau; % ms
TR = [2000; 2000; 2000; 2000; 2000]; % ms
tau_rf = 10;
flipAngle(1) = pi/2
flipAngle(2) = pi
total_time = sum(TR);

% M = zeros(3, total_time);

M0 = 1
S0 = 0;
offset = 0.000000001*3*43 %frequency offset (MHz)
time=1;
%%
for n = 1:numel(TR)
    
    
    for t =  time
        
        % M(:,t) = Rot_x(pi/2)*M0
        
        M(:,t) = M0*[0 ; sin(flipAngle(1)); cos(flipAngle(1))]
        
        Mtransverse(t) = complex(M(1,t),M(2,t));
        
       
    end
    %
    for t = time + 1 :(time + tau(n))
        
        
        % M(2,t) = M(2,1)*exp(-(t-1)/T2) ;
        % M(3,t) = M0(3) + (M(3,1) - M0(3))*(exp(-t/T1))  ;
        % S(t) = (1 - exp(-((t-1)/T1)))*exp(-(TE/T2)) ;
        %
        
        M(:,t) = M0*[exp(-(t-time)/T2)*sin(flipAngle(1))*sin(offset*(t-time));
            exp(-(t-time)/T2)*sin(flipAngle(1))*cos(offset*(t-time));
            exp(-(t-time)/T1)*cos(flipAngle(1)) + (1 - exp(-(t-time)/T1))];
        
        Mtransverse(t) = complex(M(1,t) , M(2,t));
        
       
    end
    time = t
    
    %
    for t = time
        R = [cos(flipAngle(2)), 0, -sin(flipAngle(2));
            0, 1, 0;
            sin(flipAngle(2)), 0, cos(flipAngle(2))];
        
        M(:,t) = R*M(:,(time)) ;
        Mtransverse(t) = complex( M(1,t), M(2,t));
        
        
    end
    
     start = time;
    
    for t = time + 1 : sum(TR(1:n))
      
        M(:,t) = [exp(-(t-start)/T2)*( M(1,start)*cos(offset*(t-start)) + M(2,start)*sin(offset*(t-start)) );
            exp(-(t-start)/T2)*( -M(1,start)*sin(offset*(t-start)) + M(2,start)*cos(offset*(t-start)) );
            M0*(1 - exp(-(t-start)/T1)) + M(3,(start))*exp(-(t-start)/T1)] ;
        
        
        Mtransverse(t) = complex(M(1,t), M(2,t))  ;
        
    end
    
    time = t
    M0 =  M(3,sum(TR(1:n)))
    
end
    
    



%%


    for t = TR(n-1) + 1
        t
        % M(:,t) = Rot_x(pi/2)*M0
        
        M(:,t) = M(3,TR(n-1))*[0 ; sin(flipAngle(1)); cos(flipAngle(1))];
        Mtransverse(t) = complex(M(1,t),M(2,t));
    end
    %
    for t = (TR(n-1)+2):(TR(n-1)+tau(n) -1)
        t
     
        
        % M(2,t) = M(2,1)*exp(-(t-1)/T2) ;
        % M(3,t) = M0(3) + (M(3,1) - M0(3))*(exp(-t/T1))  ;
        % S(t) = (1 - exp(-((t-1)/T1)))*exp(-(TE/T2)) ;
        %
        
        tStart = (TR(n-1)+2);
        Start = tStart - 1 ;
      
        M(:,t) = M(3,TR(n-1))*[exp(-(t-Start)/T2)*sin(flipAngle(1))*sin(offset*(t-Start));
            exp(-(t-Start)/T2)*sin(flipAngle(1))*cos(offset*(t-Start));
            exp(-(t-Start)/T1)*cos(flipAngle(1)) + (1 - exp(-(t-Start)/T1))];
        
        M(:,t)
       
        
        Mtransverse(t) = complex(M(1,t) , M(2,t));
        
    end
    
    %
    for t = tau(n)
        R = [cos(flipAngle(2)), 0, -sin(flipAngle(2));
            0, 1, 0;
            sin(flipAngle(2)), 0, cos(flipAngle(2))];
        
        M(:,t) = R*M(:,(tau-1)) ;
        Mtransverse(t) = complex( M(1,t), M(2,t));
    end
    
    
    
    for t = ((tau(n))+1) : TR(n)
        
        % M(2,t) = M(2,(TE/2))*exp(-(t-(TE/2))/T2) ;
        %
        % M(3,t) = M0(3) + ((M(3,(TE/2)) ) - M0(3))*(exp(-(t-(TE/2))/T1));
        % M(3,t) = M(3,(TE/2))*exp(-(t-(TE/2))/T1) + M0(3)*(1-exp(-(t-(TE/2))/T1));
        
        
        M(:,t) = [exp(-(t-tau(n))/T2)*( M(1,tau(n))*cos(offset*(t-tau(n))) + M(2,tau)*sin(offset*(t-tau)) );
            exp(-(t-tau(n))/T2)*( -M(1,tau)*sin(offset*(t-tau)) + M(2,tau)*cos(offset*(t-tau)) );
            M0*(1 - exp(-(t-tau)/T1)) + M(3,(tau))*exp(-(t-tau)/T1)] ;
        
        
        Mtransverse(t) = complex(M(1,t), M(2,t))  ;
        t;
    end
    
    
end
%%
figure
hold on
% plot(M(1,:))
% plot(M(2,:))
plot(M(3,:))
plot(imag(Mtransverse))
% plot(real(Mtransverse))
% plot(M0*exp(-(1:2000)/T2))
 
xlabel 'Time (ms)'
ylabel 'Magnetisation'
legend 'Mz' 'Mx' 'T2 decay' 'TE' 'Mx' 'Mx'
set(gca,'fontsize',20)
%%
defaultTR = 2000

TR = [defaultTR, defaultTR, defaultTR, defaultTR, defaultTR];
for t = 1:5000
    
    if 0 < t && t < TR(n)
        
        Mz(t) = M0*(1-exp(-t/T1))  ;
        
        Mxy(t)  = M0*exp(-t/T2) ;
        
    end
    
    for n = 1:5
        
        if n*TR(n) <= t && t < (n+1)*TR(n)
            t
            
            Mz(t) = M0*(1 - exp(-(t-(n*TR(n)))/T2))  ;
            
            Mxy(t)  = M0*(1 - exp(-n*TR(n)/T1)) * exp(-(t-n*TR(n))/T2) ;
            
        end
        
    end
    
end
%%
figure
subplot 121
plot(Mxy,'-*')
ylabel 'Mxy'
subplot 122
plot(Mz,'-*')
ylabel 'Mz'
%%

plot(1:5,s)
plot(abs(invfft(s)))

