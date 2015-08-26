clear all
T2 = 200
T1 = 200
R2 = 1/T2
R1 = 1/T1
Mz0 = 1
totalTime = 5000
M = zeros(4,numel(totalTime));
TE = [32, 132]
TR = [130, 230]
deltaOmega = 0
omega = 0
%
flipAngles(:,1) = [90,90];
flipAngles(:,2) = [180, 180]

flipAngles
A = [ -R2 deltaOmega 0 0 ...
    ; -deltaOmega -R2 omega 0 ...
    ; 0 -omega -R1 R1*Mz0 ...
    ;0 0 0 0]

M0 = [0 0 1 1]'
M(:,1) = M0

elapsedTime = 1
for n = 1:2
% 90 degree rf pulse    
M(1:3,elapsedTime)= Rot_x(flipAngles(n,1),M(1:3,elapsedTime));

% magnetisation evoution
M0 = M(:,elapsedTime)
tau = TE(n)/2;
for t = elapsedTime+1:elapsedTime+tau
M(:,t)= expm(A*(t-elapsedTime))*M0;
end
elapsedTime = t

%180 degree rf pulse
M(1:3,elapsedTime) = Rot_x(flipAngles(n,2),M(1:3,t));

% % magnetisation evolution
M0 = M(:,elapsedTime)
for t = (elapsedTime+1):(elapsedTime + TR(n)) +500
M(:,t)= expm(A*(t-elapsedTime))*M0;
end
elapsedTime = t;
end
%% Plot magnetisation ev

figure; 
hold on 
plot(M(1,:))
plot(abs(M(2,:)))
plot(M(3,:))