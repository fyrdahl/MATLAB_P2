<<<<<<< HEAD
function [Msig, M] = sim_SE_Bernstein(T1, T2, minTR, minTE, offsets,nRuns,df)
%% SIM_SE_BERNSTEIN Simulate signal for spin echo MRI sequence, allowing for variable TRs, TEs and flip angles.
%
% Function to simulate the signal expect for a "spin-echo" sequence, with  given flip angles, echo times (TEs) and
% repetition times (TRs).
% The magnetisation evolution equations in this function are from chapter on refocusing pulses (section 3.3) in 'Handbook
% of MRI pulse sequences' by Bernstein.
%
% [MSIG, M] = SIM_SE_BERNSTEIN(T1, T2, MINTR, MINTE, OFFSETS, NRUNS, DF)
%
%   T1 is a scalar. Longitudinal relaxation constant.
%
%   T2 is a scalar. Transverse relaxation constant.
%
%   MINTR is a scalar. The minimum TR hard-coded in the MRI sequence code.
%   Set to zero if OFFSETS allows for this.
%
%   MINTE is a scalar. The minimum TE hard-coded in the MRI sequence code.
%   Set to zero if OFFSETS allows for this.
%
%   OFFSETS is a 4-column matrix. Column 1: TR values, Column 2: TE values,
%   Column 3: Flip angle 1 rotation (degrees), Column 4: Flip angle 2
%   rotation (degreese).
%
%   NRUNS is a scalar. Defines how many times the code cycles through the
%   offset list.
%
%   DF is a scalar. Defines the resonance frequence offset (Hz).
%
%
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% FMRIB, Oxford, UK.
%
%
%
%% Settings and Preparation

% How many timepoints?
offsets = repmat(offsets,nRuns,1);

Msig = zeros(size(offsets,1),1);

Mzeq = 1; %total magnetisation
M = [ 0 0 Mzeq]'; % all in z-direction

%% Signal evolution for one TR
for n = 1:size(offsets,1)  
    
    %% Evolution from pulse 1 to TE/2.
    % M vector precesses about the z axis.
    % Assume Mx and My are zero. Just rotate current Mz.
    % TE has been verified for 1 slice, 'fix slice burst flag' = ON
    TE = offsets(n,2) + minTE;
    tau = TE/2;
    M =  M(3)*[exp(-tau/T2)*sind(offsets(n,3))*sind(df*tau); ...
        exp(-tau/T2)*sind(offsets(n,3))*cosd(df*tau); ...
        exp(-tau/T1)*cosd(offsets(n,3)) + (Mzeq/M(3))*(1-exp(-tau/T1))];
    % For example, if df = 0:
    %       Mx = 0
    %       My = Mz(0)*exp(-tau/T2)*sin(FA1)
    %       Mz = Mz(0)*exp(-tau/T1)*cos(FA1) + Mzeq*(1 - exp(-tau/T1))
    
    %% Effect of pulse 2 at TE/2
    %counter-clockwise rotation about the y axis, cause by pulse 2
    M = rot_y(offsets(n,4))*M;

    %% Magnetisation at TE
    Mte = [ exp(-tau/T2)*(M(1)*cosd(df*tau) + M(2)*sind(df*tau)); ...
        exp(-tau/T2)*(-M(1)*sind(df*tau) + M(2)*cosd(df*tau)); ...
        Mzeq*(1-exp(-tau/T1)) + M(3)*exp(-tau/T1)];
    % Absolute Signal at TE
    Msig(n,1) = complex(Mte(1),Mte(2));
    
    %% Magnetisation at end of TR
    % TR = TRoffset + protocolCardTR + TEoffset
    % (TR has been verified for 1 slice, 'fix slice burst flag' = ON
    TR = offsets(n,1) + minTR;
    t = TR - TE;
    M = [ exp(-t/T2)*(M(1)*cosd(df*t) + M(2)*sind(df*t)); ...
        exp(-t/T2)*(-M(1)*sind(df*t) + M(2)*cosd(df*t)); ...
        Mzeq*(1-exp(-t/T1)) + M(3)*exp(-t/T1)];
    
    
    M(1:2) = 0; % Assume the transverse signal is spoilt
    
end

=======
function [M, Mxy,flipAngles,imageTimes, t0s] = sim_SE_bernstein(T1, T2, fingerprintOffsetList,freqOffset, nSlices, nTimeCoursePts,plotFlag)
%% SimSE_Bernstein
% Function to simulate the signal expect for given flip angles, echo times (TEs) and
% repetition times (TRs).
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
%
% The magnetisation evolution equations in this function are from 'Handbook of MRI pulse sequences - section 3.3'
%%
nOffsets = size(fingerprintOffsetList,1);
if rem(nTimeCoursePts,nOffsets) ~= 0
    error('number of time points is not a multiple of the number of timings list entries')
end

nRepeats = nTimeCoursePts/nOffsets;

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
    tau = round(t0 + 1 + (TEmin + TEoffsets(n))/2) ;
    for t = [(t0 + 1), tau]
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

if strcmp(plotFlag,'plot')
    plot_simulated_signal(M, Mxy, imageTimes);
    dateTime = datestr(now,30);
    disp(['Plotted simulated signal. Date and time:',dateTime(1,1:13)]) %yyyy/mm/dd T hh:mm
end
>>>>>>> 47303176e9617d0556d5fc1882e1cda83864dfe1
end

