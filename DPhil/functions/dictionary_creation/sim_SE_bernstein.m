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

end

