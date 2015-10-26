function [M, Mxy,flipAngles, t0s] = SimBloch(T1, T2, fingerprintOffsetList, plotFlag,freqOffset, nSlices)
% Jack Allen

% see Murase 2011

clear simImageMtransverse
clear M

TRmin = 130;
TEmin = 32;
TRoffsets = fingerprintOffsetList(:,1);
TEoffsets = fingerprintOffsetList(:,2);
flipAngles(:,1) = degtorad(fingerprintOffsetList(:,3));
flipAngles(:,2) = degtorad(fingerprintOffsetList(:,4));

R2 = 1/T2;
R1 = 1/T1;
% initial magnetisation and time
Mzeq = 1; %equilibrium magnetisation
Mxy = zeros(1,numel(fingerprintOffsetList(:,1)));
t0 = 1;
M(:,t0) = [0;0;1;1]'; %initial magnetisation vector is along z

A = [-R2 0 0 0;
    0 -R2 0 0;
    0 0 -R1 R1*Mzeq;
    0 0 0 0];

%% Signal evolution

for n = 1:numel(fingerprintOffsetList(:,1))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Immediately before Pulse 1
    %     Mtransverse(t0) = complex( M(1,t0), M(2,t0));
    t0s(n) = t0;
    
    % Effect of Pulse 1
    % disp('flip 1')
    %  M(:,t0+1) = M(3,t0)*[0, sin(flipAngles(n,1)), cos(flipAngles(n,1)), 1]';
     Rot1 = [cos(flipAngles(n,1)), 0, sin(flipAngles(n,1)), 0;
        0, 1, 0, 0;
        -sin(flipAngles(n,1)), 0, cos(flipAngles(n,1)), 0;
        0, 0, 0,1];
   % Rot1 = [1, 0, 0, 0;
     %   0, cos(flipAngles(n,1)), sin(flipAngles(n,1)), 0;
     %   0, -sin(flipAngles(n,1)),cos(flipAngles(n,1)), 0;
      %  0, 0, 0,1];
      
    %     Mtransverse(t0+1) = complex( M(1,t0+1), M(2,t0+1));
    M(1:2,t0) = 0;
    M(:,t0+1) = Rot1*M(:,t0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evolution from after pulse 1, until pulse 2
    % disp('decay 1')
    
    for t = (t0 + 1) : (t0 + 1 + (TEmin + TEoffsets(n))/2)
        M(:,t) = expm(A*(t - (t0 + 1)))*M(:,t0+1);
        
        % e.g:  M(1,t) = Mx(t) = 0 if freqOffset = 0
        %       M(2,t) = My(t) = Mz(t0)*sin(FA1) if freqOffset = 0
        %       M(3,t) = Mz(t) = Mz(t0)*cos(FA1) + Mzeq*(1 - exp(-(t-elapsedTime)/T1))
        
        %         Mtransverse(t) = complex(M(1,t) , M(2,t));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Effect of pulse 2 at TE/2
    %  disp('flip 2');
    
    %rotation about the y axis
    %Rot2 = [cos(flipAngles(n,2)), 0, -sin(flipAngles(n,2)), 0;
    %    0, 1, 0, 0;
    %    sin(flipAngles(n,2)), 0, cos(flipAngles(n,2)), 0;
     %   0, 0, 0,1];
    Rot2 = [1, 0, 0, 0;
     0, cos(flipAngles(n,2)), -sin(flipAngles(n,2)), 0;
     0, sin(flipAngles(n,2)),cos(flipAngles(n,2)), 0;
      0, 0, 0,1];
    
    tau = (t0 + 1 + (TEmin + TEoffsets(n))/2) ;
    
    M(:,tau+1) = Rot2*M(:,tau);
    %     Mtransverse(tau) = complex( M(1,tau), M(2,tau));
    %     Mtransverse(tau+1) = complex( M(1,tau+1), M(2,tau+1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evolution from after pulse 2, until the newly calculated TR
    % disp('decay 2')
    
    for t = (tau+1) : t0 + TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n)
        
        M(:,t) = expm(A*(t-(tau+1)))*M(:,tau+1);
        
        
        % sample the magnetisation at TE(n)
        if t == (t0 + TEmin + TEoffsets(n))
            if Mxy(n) == 0
                Mxy(n) = abs(complex(M(1,t), M(2,t)));
                %Mxy(n) =  M(2,t);
                %             Mxy(n) = M(3,t0)*( sin(flipAngles(n,1))*sin((flipAngles(n,2))/2)*sin((flipAngles(n,2))/2)*exp(-(t - t0)/T2) );
                imageTimes(n) = t;
                t0s(n) = t0;
            end
        end
  
    end
    
    % newly calculated newTRs
    newTRs(n) = TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n);
    
    % update total time that has passed
    t0 = t0 + newTRs(n)
    
    
end
%% plotting the signal evolution
if plotFlag == 'showPlot'
    
    figure
    hold on
    plot(M(1,:))
    plot(M(2,:))
    plot(M(3,:))
    plot(imageTimes, Mxy,'+')
    xlabel 'Time (ms)'
    ylabel 'Magnetisation'
    legend 'Mx' 'My' 'Mz' 'Mtransverse'
    set(gca,'fontsize',20)
    
end

% disp('bloch simulation: complete')
end
