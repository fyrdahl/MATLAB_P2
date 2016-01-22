function [M, Mxy,imageTimes,flipAngles, t0s] = SimSE(T1, T2, fingerprintOffsetList,nRepeats,df, nSlices, plotFlag)
%% Jack Allen 
% Uses solutions to the Bloch equations to simulate the expected signal for
% a spin echo experiment. The tutorial by Brian Hargreaves was used to
% guide the writing of this code
% (http://www-mrsrl.stanford.edu/~brian/bloch/).

%%
TRmin = 130;
TEmin = 32;
TRoffsets = fingerprintOffsetList(:,1);
TEoffsets = fingerprintOffsetList(:,2);
flipAngles(:,1) = fingerprintOffsetList(:,3);
flipAngles(:,2) = fingerprintOffsetList(:,4);
nTimePts = nRepeats*numel(fingerprintOffsetList(:,1));
nImages = numel(fingerprintOffsetList(:,1));

imageTimes = zeros(nTimePts,1);
newTRs = zeros(nTimePts,1); %
t0s = zeros(nTimePts,1); %To track the total time elapsed
Mzeq = 1; %equilibrium magnetisation
Mxy = zeros(1,nTimePts);
t0 = 1; %initial time point
M(:,t0) = [0;0;Mzeq]'; %initial magnetisation vector is along z
dt = 1; %time step size

%function from Stanford tutorial website to simulate magnetisation evolution
[A, B] =  freeprecess(dt, T1, T2, df);

%% Signal evolution
for i = 1:nRepeats
    for n = 1:nImages

      
        imageIndex = n + nImages*(i - 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Effect of Pulse 1
        Rot_y(flipAngles(n,1));
        M(:,t0) = Rot_y(flipAngles(n,1))*M(:,t0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % evolution from after pulse 1, until pulse 2
        % disp('decay 1')
        for t = (t0 + 1) : (t0 + 1 + (TEmin + TEoffsets(n))/2)
            M(:,t) = A*M(:,t-1) + B;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Effect of pulse 2 at TE/2
        %  disp('flip 2');
        tau = (t0 + 1 + (TEmin + TEoffsets(n))/2) ;
        %rotation about the x axis
        M(:,tau) = Rot_x(flipAngles(n,2),M(:,tau));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % evolution from after pulse 2, until the newly calculated TR
        % disp('decay 2')
        for t = (tau+1) : t0 + TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n)
            M(:,t) = A*M(:,t-1) + B;
            % sample the magnetisation at TE(n)
            if t == (t0 + TEmin + TEoffsets(n))             
                if Mxy(imageIndex) == 0
                 imageTimes(imageIndex) = t;
                  Mxy(1,imageIndex) = abs(complex( M(1,t), M(2,t) ));    
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % newly calculated newTRs
        newTRs(imageIndex) = TRmin + (nSlices*TEoffsets(n)) + TRoffsets(n);
        
        % update total time that has passed
        t0 = t0 + newTRs(n);
        M(1:2,t0) = 0;
        
        t0s(imageIndex) = t0;
    end
    
%% plotting the signal evolution
if plotFlag == 1
figure
hold on
plot(M(1,:))
plot(M(2,:))
plot(M(3,:))
plot(abs(complex(M(1,:),M(2,:))),'--k')
plot(imageTimes, Mxy,'+')
xlabel 'Time (ms)'
ylabel 'Magnetisation'
legend 'Mx' 'My' 'Mz' 'Mtransverse'
set(gca,'fontsize',20)  
end

end
