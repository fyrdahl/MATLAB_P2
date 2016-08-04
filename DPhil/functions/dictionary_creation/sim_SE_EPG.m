function [sig, P] = sim_SE_EPG(RFpulses1,RFpulses2, TRs, TEs, T1, T2)
% SIM_SE_EPG counter-clockwise rotation (degrees) about the x-axis
%
% Function to simulate (using Extended Phase Graphs EPGs) the signal from a spin-echo EPI sequence. All times must in seconds and angles in
% Radians. Code adapted from FISP code from Lior Weizman (which was adapted
% from Brian Hargeaves' (Stanford) code: http://web.stanford.edu/class/rad229/Matlab.html. )
%
% [SIG] = SIM_SE_EPG(RFPULSES1, RFPULSES2, TRS, TE, T1, T2)
%
%   RFPULSES1 is a complex vector. Represents all of the excitation pulses
%   Absolute values are rotations in radians.
%
%   RFPULSES2 is a complex vector. Represents all of the 'refocusing'
%   pulses. Absolute values are rotations in radians.
%
%   TRS is a vector. Contains the length of each TR (in seconds)
%
%   TES is a vector. Contains the length of each TE (in seconds).
%   
%   T1 is a scalar. Defines the longitudinal relaxation constant (in
%   seconds).
%
%   T2 is a scale. Defines the transverse relaxation constant (in
%   seconds).
%
%   SIG is a complex vector. Contains the simulated signal for each readout (i.e.
%   each TR).
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% Copyright © 2016 University of Oxford

Ntr = length(TRs);	% Number of TRs
Nstates = 20;	% Number of states to simulate 

P = zeros(3,Nstates);	% State matrix
P(3,1)=1;		% Equilibrium magnetization.

sig=zeros(1,length(TRs)); % Vector holding the recieved signals

FA1=abs(RFpulses1);
FA2=abs(RFpulses2);
FP1=angle(RFpulses1);
FP2=angle(RFpulses2);

for n = 1:Ntr
   
    % RF Excitation pulse
    P = epg_rf(P,FA1(n),FP1(n));
      
    % Relxation from RF until 'refocusing' pulse at TE/2
    P = epg_grelax(P,T1,T2,TEs(n)/2,0,0,0,0);
    
    % RF 'refocusing' pulse
    P = epg_rf(P,FA2(n),FP2(n));
    
    % Relaxation from time TE/2 to TE
    P = epg_grelax(P,T1,T2,TEs(n)/2,0,0,0,0); 

    % Signal is F0 state.
    sig(n) = P(1,1);
     
    % relaxation and spoiler gradient
    P = epg_grelax(P,T1,T2,TRs(n)-TEs(n),1,0,1,0);
    
end;

end
