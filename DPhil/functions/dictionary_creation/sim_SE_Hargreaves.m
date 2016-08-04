function [Msig, M] = sim_SE_Hargreaves(T1, T2, minTR, minTE, offsets,nRuns, df, plotFlag)
% SIM_SE Spin Echo (SE) Bloch simulations
%
% Function that uses solutions to the Bloch equations to simulate the expected signal for
% a spin echo experiment. The tutorial by Brian Hargreaves was used to
% guide the writing of this code
% (http://www-mrsrl.stanford.edu/~brian/bloch/).
%
% [M, MXY,IMAGETIMES,FLIPANGLES, T0S] = SIM_SE(T1, T2, FINGERPRINTOFFSETLIST, NREPEATS, DF, NSLICES, PLOTFLAG)
%
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
%
%




nImages = nRuns*size(offsets,1);
Msig = zeros(1,nImages);
M = zeros(3,nImages);
Mzeq = 1;
M0 = [ 0; 0; Mzeq];
offsets = repmat(offsets,[nRuns,1]);

%% Signal evolution

for n = 1:nImages
    
    
    % Flip 1, a rotation about y-axis
    M = rot_y(offsets(n,3))*M0;
    % Magnetisation evolution over time until TE/2
    %(TE checked for 1 slice, 'fix slice burst flag' = ON)
    TE = offsets(n,2) + minTE;
    
    %evolution from FA1 to TE/2
    % calculate relaxation coefficients (taken from Hargeaves 'freeprecess' function
    fp = [exp(-TE/2/T2) 0 0; ...
        0 exp(-TE/2/T2) 0; ...
        0 0 exp(-TE/2/T1)]*zrot(2*pi*df*TE/2/1000); %rotate by Resonant precession, radians.
    M(1) = M(1)*fp(1,1);
    M(2) = M(2)*fp(2,2);
    M(3) = Mzeq - (Mzeq - M(3))*fp(3,3);
    
    % Flip 2, a rotation about the x axis
    M = rot_y(offsets(n,4))*M;
    
    % Magnetisation evolution over time until TE
    M(1) = M(1)*fp(1,1);
    M(2) = M(2)*fp(2,2);
    M(3) = Mzeq - (Mzeq - M(3))*fp(3,3);
    
    % Collect Signal
    Msig(n) = complex(M(1),M(2));
    
    %Evolution from TE to TR
    %TE to TR = minTR + TRadd + TEadd - (TEadd + minTE)
    %(TR checked for 1 slice, 'fix slice burst flag' = ON)
    t = minTR+offsets(n,1)-minTE - TE;
    fp = [exp(-t/T2) 0 0;0 exp(-t/T2) 0;0 0 exp(-t/T1)]*zrot(2*pi*df*t/1000); %rotate by Resonant precession, radians.
    M(1) = M(1)*fp(1,1);
    M(2) = M(2)*fp(2,2);
    M(3) = Mzeq - (Mzeq - M(3))*fp(3,3);
    
    % At the end of the TR, set the current Mz to be the initial Mz for the next TR
    M0 = [ 0;  0; M(3)]; % ignore the remaining transverse signal
end

%% plotting the signal evolution
if strcmp(plotFlag,'plot')
    figure
    hold on
    plot(M(1,:))
    plot(M(2,:))
    plot(M(3,:))
    plot(Msig,'+')
    xlabel 'Time (ms)'
    ylabel 'Magnetisation'
    legend 'Mx' 'My' 'Mz' 'Mtransverse'
    set(gca,'fontsize',20)
end

end
