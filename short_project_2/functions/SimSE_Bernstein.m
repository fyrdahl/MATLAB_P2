function [M, Mxy,flipAngles,imageTimes, t0s] = SimSE_Bernstein(T1, T2, fingerprintOffsetList, plotFlag,freqOffset, nSlices, nTimeCoursePts)
% Jack Allen

% magnetisation evolution equations in this function are from 'Handbook of MRI pulse sequences - section 3.3'
%%

if rem(nTimeCoursePts,24) ~= 0
    error('number of time points is not a multiple of the number of timings list entries')
end

nRepeats = nTimeCoursePts/24;

TRmin = 130;
TEmin = 32;
TRoffsets = fingerprintOffsetList(:,1);
TEoffsets = fingerprintOffsetList(:,2);
TRoffsets = repmat(TRoffsets,[nRepeats,1]);
TEoffsets = repmat(TEoffsets,[nRepeats,1]);

flipAngles(:,1) = degtorad(fingerprintOffsetList(:,3));
flipAngles(:,2) = degtorad(fingerprintOffsetList(:,4));
flipAngles = repmat(flipAngles,[nRepeats,1]);


Mzeq = 1; %equilibrium magnetisation
t0 = 1;
M(:,t0) = [0;0;1]'; %initial magnetisation vector is along z
Mxy = zeros(1, nTimeCoursePts);

%% Signal evolution
for n = 1:nTimeCoursePts
    
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
    
    M(:,tau) = R*M(:,tau);
    %     Mtransverse(tau) = complex( M(1,tau), M(2,tau));
    %     Mtransverse(tau+1) = complex( M(1,tau+1), M(2,tau+1));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evolution from after pulse 2, until the newly calculated TR
    % disp('decay 2')
    for t = (tau) : t0 + TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n)
        
        M(:,t) = [exp(-(t-tau)/T2)*( M(1,tau)*cos(freqOffset*(t-tau)) + M(2,tau)*sin(freqOffset*(t-tau)) );
            exp(-(t-tau)/T2)*( -M(1,tau)*sin(freqOffset*(t-tau)) + M(2,tau)*cos(freqOffset*(t-tau)) );
            Mzeq*(1 - exp(-(t-tau)/T1)) + M(3,tau)*exp(-(t-tau)/T1)] ;
        
        %         Mtransverse(t) = complex(M(1,t), M(2,t))  ;
        
        % sample the magnetisation at TE(n)
        if t == (t0 + TEmin + TEoffsets(n))
            
            %Mxy(n) = abs(Mtransverse(t));
            Mxy(n) = abs(M(3,t0)*( sin(flipAngles(n,1))*sin((flipAngles(n,2))/2)*sin((flipAngles(n,2))/2)*exp(-(t - t0)/T2) ));
            
            imageTimes(n) = t;
            t0s(n) = t0;
        end
    end
    
    % newly calculated newTRs
    newTRs(n) = TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n);
    
    % update total time that has passed
    t0 = t0 + newTRs(n);
    
end
end

