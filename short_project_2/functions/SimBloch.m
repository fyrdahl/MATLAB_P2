function [M, Mxy,flipAngles, t0s] = SimBloch(T1, T2, fingerprintOffsetList, plotFlag,freqOffset, nSlices, nTimeCoursePts)
% Jack Allen

% magnetisation evolution equations in this function are from 'Handbook of MRI pulse sequences - section 3.3'

clear simImageMtransverse
clear M

TRmin = 130;
TEmin = 32;
TRoffsets = fingerprintOffsetList(:,1);
TEoffsets = fingerprintOffsetList(:,2);
flipAngles(:,1) = degtorad(fingerprintOffsetList(:,3));
flipAngles(:,2) = degtorad(fingerprintOffsetList(:,4));

% initial magnetisation and time
Mzeq = 1; %equilibrium magnetisation
t0 = 1;
M(:,t0) = [0;0;1]'; %initial magnetisation vector is along z
Mxy = zeros(1, numel(fingerprintOffsetList(:,1)));

%% Signal evolution

for n = 1:(numel(fingerprintOffsetList(:,1)))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Immediately before Pulse 1
%     Mtransverse(t0) = complex( M(1,t0), M(2,t0));
    t0s(n) = t0;
    
    % Effect of Pulse 1
    % disp('flip 1')
    M(:,t0+1) = M(3,t0)*[0, sin(flipAngles(n,1)), cos(flipAngles(n,1))]';
%     Mtransverse(t0+1) = complex( M(1,t0+1), M(2,t0+1));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evolution from after pulse 1, until pulse 2
    % disp('decay 1') 
  
     for t = (t0 + 1) : (t0 + 1 + (TEmin + TEoffsets(n))/2)
        M(:,t) =  M(3,t0)*[exp(-(t-t0 + 1)/T2)*sin(flipAngles(n,1))*sin(freqOffset*(t-t0 + 1));
            exp(-(t-t0 + 1)/T2)*sin(flipAngles(n,1))*cos(freqOffset*(t-t0 + 1));
            exp(-(t-t0 + 1)/T1)*cos(flipAngles(n,1)) + (Mzeq/M(3,t0))*(1 - exp(-(t-t0 + 1)/T1))];
  
        % e.g:  M(1,t) = Mx(t) = 0 if freqOffset = 0
        %       M(2,t) = My(t) = Mz(t0)*sin(FA1) if freqOffset = 0
        %       M(3,t) = Mz(t) = Mz(t0)*cos(FA1) + Mzeq*(1 - exp(-(t-elapsedTime)/T1))
        
%         Mtransverse(t) = complex(M(1,t) , M(2,t));
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Effect of pulse 2 at TE/2
    %  disp('flip 2');
    
    %rotation about the y axis
    R = [cos(flipAngles(n,2)), 0, -sin(flipAngles(n,2));
        0, 1, 0;
        sin(flipAngles(n,2)), 0, cos(flipAngles(n,2))];
    
    tau = (t0 + 1 + (TEmin + TEoffsets(n))/2) ;
    
    M(:,tau+1) = R*M(:,tau);
%     Mtransverse(tau) = complex( M(1,tau), M(2,tau));
%     Mtransverse(tau+1) = complex( M(1,tau+1), M(2,tau+1));
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evolution from after pulse 2, until the newly calculated TR
    % disp('decay 2')
    
   for t = (tau+1) : t0 + TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n)
   
        M(:,t) = [exp(-(t-tau + 1)/T2)*( M(1,tau + 1)*cos(freqOffset*(t-tau)) + M(2,tau + 1)*sin(freqOffset*(t-tau)) );
            exp(-(t-tau + 1)/T2)*( -M(1,tau + 1)*sin(freqOffset*(t-tau + 1)) + M(2,tau + 1)*cos(freqOffset*(t-tau)) );
            Mzeq*(1 - exp(-(t-tau + 1)/T1)) + M(3,tau + 1)*exp(-(t-tau + 1)/T1)] ;
        
             
        
%         Mtransverse(t) = complex(M(1,t), M(2,t))  ;
        
        % sample the magnetisation at TE(n)
        if t == (t0 + TEmin + TEoffsets(n))
            if Mxy(n) == 0
                n;
%                 Mxy(n) = abs(Mtransverse(t));

   
                Mxy(n) = abs(M(3,t0)*( sin(flipAngles(n,1))*sin((flipAngles(n,2))/2)*sin((flipAngles(n,2))/2)*exp(-(t - t0)/T2) ));
   
               % if Mxy(n) < 0
                %    disp('NEGATIVE')
                 %   pause
               % end
                
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
