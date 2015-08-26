function [Mtransverse, simImageMtransverse, M,Mxy,flipAngles] = SimBloch(T1, T2, fingerprintOffsetList, plotFlag,freqOffset)

%  magnetisation evolution equations from 'Handbook of MRI pulse sequences - section 3.3'
%
clear simImageMtransverse
clear M

deltaOmega = 0;
omega = 0;

TRmin = 130;
TEmin = 32;
TRoffsets = fingerprintOffsetList(:,1);
TEoffsets = fingerprintOffsetList(:,2);

flipAngles(:,1) = degtorad(fingerprintOffsetList(:,3));
flipAngles(:,2) = degtorad(fingerprintOffsetList(:,4));

nSlices = 2;

t0 = 1;
Mzeq = 1;
M0 = 1;
elapsedTime = 1;
M(:,t0) = M0*[0;0;1]';

%%

for n = 1:numel(fingerprintOffsetList(:,1))
    
    t0 = elapsedTime;
    t0s(n) = t0
    %     disp('flip 1');
    
    %slice 1
    M(:,t0+1) = M0*[0, sin(flipAngles(n,1)), cos(flipAngles(n,1))]';
    
    Mtransverse(t0) = complex( M(1,t0), M(2,t0));
    Mtransverse(t0+1) = complex( M(1,t0+1), M(2,t0+1));
    
    % evolution after first flip pulse
    % disp('decay 1');
    elapsedTime = t0 + 1;
    for t = (elapsedTime) : (elapsedTime + (TEmin + TEoffsets(n))/2)
        
        % M(1,:) = 0 if freqOffset = 0       
        % M(2,:) = M(3,t0)*sin(FA1) if freqOffset = 0
        % M(3,:) = M(3,t0)*cos(FA1) + Mzeq*(1 - exp(-(t-elapsedTime)/T1))
        M(:,t) =  M(3,t0)*[exp(-(t-elapsedTime)/T2)*sin(flipAngles(n,1))*sin(freqOffset*(t-elapsedTime)); 
            exp(-(t-elapsedTime)/T2)*sin(flipAngles(n,1))*cos(freqOffset*(t-elapsedTime));
            exp(-(t-elapsedTime)/T1)*cos(flipAngles(n,1)) + (Mzeq/M(3,t0))*(1 - exp(-(t-elapsedTime)/T1))];
        
        Mtransverse(t) = complex(M(1,t) , M(2,t));
    end
    elapsedTime = t;
    % TE/2  second flip
    %  disp('flip 2');
    
    %rotation about the x axis
    %         R = [1, 0, 0;
    %             0, cos(flipAngles(n,2)), -sin(flipAngles(n,2));
    %             0, sin(flipAngles(n,2)), cos(flipAngles(n,2))];
    
    %rotation about the y axis
    R = [cos(flipAngles(n,2)), 0, -sin(flipAngles(n,2));
        0, 1, 0;
        sin(flipAngles(n,2)), 0, cos(flipAngles(n,2))];
        
    M(:,elapsedTime+1) = R*M(:,elapsedTime);
    Mtransverse(elapsedTime) = complex( M(1,elapsedTime), M(2,elapsedTime));
    
    tau = elapsedTime + 1;
    
    % evolution from TE/2 to TR
    % disp('decay 2')
    for t = (tau) : t0 + TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n)
        
        M(:,t) = [exp(-(t-tau)/T2)*( M(1,tau)*cos(freqOffset*(t-tau)) + M(2,tau)*sin(freqOffset*(t-tau)) );
            exp(-(t-tau)/T2)*( -M(1,tau)*sin(freqOffset*(t-tau)) + M(2,tau)*cos(freqOffset*(t-tau)) );
            Mzeq*(1 - exp(-(t-tau)/T1)) + M(3,tau)*exp(-(t-tau)/T1)] ;
                        
        Mtransverse(t) = complex(M(1,t), M(2,t))  ;
        
        if t == (t0 + TEmin + TEoffsets(n))
            n;
            simImageMtransverse(n) = abs(Mtransverse(t));
            Mxy(n) = M(3,t0)*( sin(flipAngles(n,1))*sin((flipAngles(n,2))/2)*sin((flipAngles(n,2))/2)*exp(-(t-t0)/T2) );
            imageTimes(n) = t
            t0s(n) = t0;
        end
        
    end
    newTRs(n) = TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n)
    %     newTRs(n) = (TRmin/2) + nSlices*TEoffsets(n) + TEmin + TRoffsets(n)
    
    elapsedTime = t0 + newTRs(n);
    
    % Mz at TR, which will be used as M0 for the next loop
    M0 =  M(3,elapsedTime);
    
end

if plotFlag == 'showPlot'
    %% plot signal evolution
    
    figure
    hold on
    plot(M(1,:))
    plot(M(2,:))
    plot(M(3,:))
%     plot(abs(Mtransverse) ,'--')
%     plot(imageTimes, simImageMtransverse,'*')
    plot(imageTimes, Mxy,'+')
    xlabel 'Time (ms)'
    ylabel 'Magnetisation'
    legend 'Mx' 'My' 'Mz' 'absolute Mtransverse'
    set(gca,'fontsize',20)
    
end

%%%%%% % % %%%
disp('bloch simulation: complete')
end
